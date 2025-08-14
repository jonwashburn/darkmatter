#!/usr/bin/env python3
"""
Benchmark (global-only, no per-galaxy tuning): RS kernel vs MOND

- Reuses SPARC master table (results/sparc_master.pkl)
- Uses the same global sigma model/masks as ilg_pure_solver
- RS: accel kernel with locked defaults from ilg_pure_solver
- MOND: global-only, simple mu function; single global disk M/L (1.0)

Outputs:
  results/bench_global_summary.csv
  results/bench_rs_per_galaxy.csv
  results/bench_mond_per_galaxy.csv
"""
from __future__ import annotations
import importlib.util
import pickle
from pathlib import Path
from typing import Dict, Any, Tuple
import numpy as np
import csv

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'
RESULTS.mkdir(parents=True, exist_ok=True)

def load_ilg():
    mod_path = ROOT / 'scripts' / 'ilg_pure_solver.py'
    spec = importlib.util.spec_from_file_location('ilg', str(mod_path))
    ilg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ilg)  # type: ignore
    return ilg


def mond_nu(y: np.ndarray) -> np.ndarray:
    """Simple MOND interpolation via nu-function: nu(y) = 0.5 + sqrt(0.25 + 1/y)."""
    y = np.asarray(y, float)
    return 0.5 + np.sqrt(0.25 + np.where(y > 0, 1.0 / y, 0.0))


def mond_v_model_kms(ilg, vgas: np.ndarray, vdisk: np.ndarray, vbul: np.ndarray, r_kpc: np.ndarray, ml_disk: float) -> np.ndarray:
    """Compute MOND circular velocity [km/s] from baryon components using simple nu.
    Uses a0 from ilg (A0_RS)."""
    # Effective baryon rotation from components
    vdisk_scaled = np.asarray(vdisk, float) * np.sqrt(max(ml_disk, 0.0))
    v2_b = np.maximum(vgas, 0.0) ** 2 + np.maximum(vdisk_scaled, 0.0) ** 2 + np.maximum(vbul, 0.0) ** 2
    v_b = np.sqrt(np.maximum(v2_b, 0.0))  # km/s
    # Newtonian g
    v_ms = v_b * ilg.KMS_TO_MS
    r_m = np.asarray(r_kpc, float) * ilg.KPC_TO_M
    gN = np.divide(v_ms * v_ms, np.maximum(r_m, 1.0), out=np.zeros_like(v_ms), where=r_m > 0)
    y = np.maximum(gN / ilg.A0_RS, 1e-12)
    g = gN * mond_nu(y)
    v = np.sqrt(g * r_m) / ilg.KMS_TO_MS
    return v


def evaluate_mond(master_table: Dict[str, Any], ilg, use_mask: bool = True, ml_disk: float = 1.0) -> Tuple[list, Dict[str, Any]]:
    """Evaluate MOND per-galaxy with the same sigma/masks as ilg; return per-galaxy rows and summary."""
    rows = []
    w_all = []
    for name, gal in master_table.items():
        df = gal.get('data')
        r_full = np.asarray(gal.get('r', []), float)
        v_obs_full = np.asarray(gal.get('v_obs', []), float)
        if r_full.size < 2 or v_obs_full.size != r_full.size:
            continue
        # beam mask
        distance_mpc = float(gal.get('distance', 10.0))
        beam_kpc = ilg.BEAM_ARCSEC * (distance_mpc * 1.0e3) / 206265.0 if use_mask else 0.0
        mask = r_full >= beam_kpc
        r = r_full[mask]
        v_obs = v_obs_full[mask]
        if r.size < 2:
            continue
        # baryon components
        if hasattr(df, 'columns') and all(c in df.columns for c in ('vgas','vdisk','vbul','verr')):
            vgas = np.asarray(df['vgas'].values, float)[mask]
            vdisk = np.asarray(df['vdisk'].values, float)[mask]
            vbul = np.asarray(df['vbul'].values, float)[mask]
        else:
            # fallback to provided v_baryon
            vb = np.asarray(gal.get('v_baryon', []), float)[mask]
            vgas = np.zeros_like(vb); vdisk = np.zeros_like(vb); vbul = vb
        # galaxy type for drift
        v_max = float(np.max(v_obs)) if v_obs.size else 0.0
        gtype = 'dwarf' if v_max < 80.0 else 'spiral'
        # sigma eff
        seff = ilg.sigma_components(df, r, v_obs, distance_mpc, gtype, float(gal.get('R_d', 3.0)), idx_mask=mask)
        # MOND velocity and chi2
        v_model = mond_v_model_kms(ilg, vgas, vdisk, vbul, r, ml_disk)
        resid = v_obs - v_model
        chi2 = float(np.sum(((resid) / seff) ** 2))
        rows.append({
            'name': name,
            'N': int(v_obs.size),
            'chi2_reduced': chi2 / max(1, int(v_obs.size)),
            'v_obs_mean': float(np.nanmean(v_obs)),
            'resid_rms': float(np.sqrt(np.nanmean(resid**2)))
        })
        # fake w_total placeholder for parity
        w_all.append(np.ones_like(v_obs))
    chis = np.array([r['chi2_reduced'] for r in rows]) if rows else np.array([np.nan])
    summary = {
        'kernel': 'MOND',
        'ml_disk': ml_disk,
        'N_galaxies': len(rows),
        'chi2_reduced_median': float(np.nanmedian(chis)),
        'chi2_reduced_mean': float(np.nanmean(chis))
    }
    return rows, summary


def main():
    ilg = load_ilg()
    master = pickle.load(open(RESULTS / 'sparc_master.pkl', 'rb'))

    # RS detailed (locked defaults)
    rs_rows, rs_summary = ilg.evaluate_detailed(master, kernel='accel', use_n_profile=True, use_xi=True, ml_disk=1.0)

    # MOND detailed (global-only)
    mond_rows, mond_summary = evaluate_mond(master, ilg, use_mask=True, ml_disk=1.0)

    # Save per-galaxy CSVs
    with open(RESULTS / 'bench_rs_per_galaxy.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['name','N','chi2_reduced','v_obs_mean','resid_rms'])
        for r in sorted(rs_rows, key=lambda x: x['chi2_reduced'], reverse=True):
            w.writerow([r['name'], r['N'], r['chi2_reduced'], r['v_obs_mean'], r['resid_rms']])
    with open(RESULTS / 'bench_mond_per_galaxy.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['name','N','chi2_reduced','v_obs_mean','resid_rms'])
        for r in sorted(mond_rows, key=lambda x: x['chi2_reduced'], reverse=True):
            w.writerow([r['name'], r['N'], r['chi2_reduced'], r['v_obs_mean'], r['resid_rms']])

    # Save summary CSV
    with open(RESULTS / 'bench_global_summary.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['model','N_gal','chi2_median','chi2_mean'])
        w.writerow(['RS_global', rs_summary['N_galaxies'], rs_summary['chi2_reduced_median'], rs_summary['chi2_reduced_mean']])
        w.writerow(['MOND_global', mond_summary['N_galaxies'], mond_summary['chi2_reduced_median'], mond_summary['chi2_reduced_mean']])

    print('Saved:', RESULTS / 'bench_global_summary.csv')
    print('Saved:', RESULTS / 'bench_rs_per_galaxy.csv')
    print('Saved:', RESULTS / 'bench_mond_per_galaxy.csv')


if __name__ == '__main__':
    main()


