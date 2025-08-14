#!/usr/bin/env python3
"""
ILG Pure Solver (no per-galaxy fitting)

- Uses globally fixed constants: alpha=0.191
- Small-lag model: w_t(r) = 1 + c_lag * [(T_dyn/T_ref)^alpha - 1]
  with c_lag = phi^{-5} from RS coherence/dual-balance scaling
- Acceleration kernel: w_g(r) = 1 + c_lag * [(g_bar/a0)^(-alpha) - 1]
  with a0 = 1.2e-10 m/s^2 (RS-derived global scale)
- Global analytic radial profile: n(r) = 1 + A*(1 - exp(-(r/r0)^p))
  with (A, r0 [kpc], p) = (7.0, 8.0, 1.6) fixed once
- n(r) is normalised so that its universal disc-weighted mean equals 1
  under W(r) ∝ r * exp(-r/R_w) with R_w=3 kpc
- Optional global complexity factor: xi = 1 + c_xi * f_gas_true^{gamma_xi}
  with c_xi = phi^{-5}, gamma_xi = 1/2, no per-galaxy tuning
- Vertical disk correction: zeta(r) with global h_z/R_d
- Chi^2 uses a standard uncertainty floor + beam smearing + asymmetric drift:
  sigma_eff^2 = verr^2 + sigma0^2 + (f*v_obs)^2 + sigma_beam^2 + sigma_asym^2
- Optional global disk M/L scale (ml_disk) applied to vdisk uniformly
- v_model^2(r) = w_total(r) * v_baryon_eff^2(r)

Reads: results/sparc_master.pkl
Writes: results/ledger_final_pure_results.pkl (bundle with variants),
        results/ilg_variants_summary.csv,
        results/ilg_best_median_per_galaxy.csv,
        results/ilg_best_mean_per_galaxy.csv
"""
from __future__ import annotations
import pickle
from pathlib import Path
from typing import Dict, Any, Tuple
import numpy as np
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Global constants
ALPHA = 0.191
PHI = (1.0 + 5.0 ** 0.5) / 2.0
C_LAG = PHI ** (-5.0)  # ~0.09017
A_N = 7.0
R0_KPC = 8.0
P_SHAPE = 1.6
R_WEIGHT_KPC = 3.0
R_MAX_KPC = 50.0
C_XI = PHI ** (-5.0)
GAMMA_XI = 0.5
A0_RS = 1.2e-10  # m/s^2
KPC_TO_M = 3.085677581491367e19
KMS_TO_MS = 1.0e3
SIGMA0_KMS = 10.0  # velocity error floor [km/s]
FRAC_FLOOR = 0.05  # fractional systematic on v_obs
BEAM_ARCSEC = 15.0
ALPHA_BEAM = 0.3  # global, no tuning
DRIFT_DWARF = 0.10  # fractional non-circular (asymmetric drift) for dwarfs
DRIFT_SPIRAL = 0.05 # fractional non-circular for spirals
HZ_RATIO = 0.25   # global h_z / R_d
ML_DISK_VALUES = [1.0, 1.2]  # test no-scale and modest global disk M/L
K_TURB = 0.07  # global outer-disk turbulence fraction (locked)
P_TURB = 1.3   # global radial shape exponent for turbulence (locked)
G_EXT_TESTS = [0.0, 0.03, 0.05]  # external field fraction of a0 (global)
ZCAP_MIN = 0.8
ZCAP_MAX = 1.2


def collect_T_ref(master_table: Dict[str, Any]) -> float:
    """Compute a global reference dynamical time (years) as median over all points."""
    all_T = []
    for gal in master_table.values():
        T = gal.get("T_dyn")
        if T is None:
            continue
        tvals = np.asarray(T, dtype=float)
        if tvals.size:
            all_T.append(tvals)
    if not all_T:
        return 1.0e8
    cat = np.concatenate(all_T)
    cat = cat[np.isfinite(cat) & (cat > 0)]
    return float(np.median(cat)) if cat.size else 1.0e8


def w_small_lag_centered(T_dyn_years: np.ndarray, T_ref_years: float) -> np.ndarray:
    ratio = np.asarray(T_dyn_years, dtype=float) / float(T_ref_years)
    time_term = np.power(np.maximum(ratio, 1e-12), ALPHA)
    return 1.0 + C_LAG * (time_term - 1.0)


def w_accel_centered(v_baryon_kms: np.ndarray, r_kpc: np.ndarray, g_ext_frac: float = 0.0) -> np.ndarray:
    v_ms = np.asarray(v_baryon_kms, dtype=float) * KMS_TO_MS
    r_m = np.asarray(r_kpc, dtype=float) * KPC_TO_M
    g_bar = np.divide(v_ms * v_ms, np.maximum(r_m, 1.0), out=np.zeros_like(v_ms), where=r_m > 0)
    g_ext = float(max(g_ext_frac, 0.0)) * A0_RS
    y = np.maximum((g_bar + g_ext) / A0_RS, 1e-12)
    y0 = np.maximum((1.0 + float(max(g_ext_frac, 0.0))), 1e-12)
    # centered so w=1 at g_bar=a0 when including external field
    return 1.0 + C_LAG * (np.power(y, -ALPHA) - np.power(y0, -ALPHA))


def n_analytic(r_kpc: np.ndarray) -> np.ndarray:
    r = np.asarray(r_kpc, dtype=float)
    base = 1.0 + A_N * (1.0 - np.exp(-np.power(np.maximum(r, 0.0) / R0_KPC, P_SHAPE)))
    return base


# Precompute global <n>_W = 1 by normalisation
_r_grid = np.linspace(0.0, R_MAX_KPC, 2001)
_r_grid[0] = 1e-6
_W = _r_grid * np.exp(-_r_grid / R_WEIGHT_KPC)
_n_vals = n_analytic(_r_grid)
_nbar = float(np.trapz(_n_vals * _W, _r_grid) / np.trapz(_W, _r_grid))
N_GLOBAL = _nbar if _nbar > 0 else 1.0


def n_normalised(r_kpc: np.ndarray) -> np.ndarray:
    return n_analytic(r_kpc) / N_GLOBAL


def xi_global(gal: Dict[str, Any]) -> float:
    f_gas = gal.get("f_gas_true")
    if f_gas is None:
        return 1.0
    try:
        val = float(f_gas)
    except Exception:
        return 1.0
    val = np.clip(val, 0.0, 1.0)
    return float(1.0 + C_XI * (val ** GAMMA_XI))


def zeta_vertical(r_kpc: np.ndarray, R_d_kpc: float) -> np.ndarray:
    r = np.asarray(r_kpc, float)
    R_d = float(R_d_kpc) if R_d_kpc and R_d_kpc > 0 else 3.0
    h_z = HZ_RATIO * R_d
    x = r / R_d
    f_profile = np.exp(-x / 2.0) * (1.0 + x / 3.0)
    zeta = 1.0 + 0.5 * (h_z / (r + 0.1 * R_d)) * f_profile
    return np.clip(zeta, ZCAP_MIN, ZCAP_MAX)


def sigma_components(df, r_kpc: np.ndarray, v_obs_kms: np.ndarray, distance_mpc: float, galaxy_type: str, R_d_kpc: float, idx_mask: np.ndarray | None = None) -> np.ndarray:
    if hasattr(df, 'columns') and 'verr' in df.columns:
        verr_full = np.asarray(df['verr'].values, float)
        verr = verr_full[idx_mask] if idx_mask is not None else verr_full
    else:
        verr = np.full_like(v_obs_kms, 5.0)
    v_obs = np.asarray(v_obs_kms, float)
    r = np.asarray(r_kpc, float)
    # base floor
    sigma_base = np.sqrt(np.maximum(verr, 0.0) ** 2 + (SIGMA0_KMS ** 2) + (FRAC_FLOOR * np.maximum(v_obs, 0.0)) ** 2)
    # beam smearing
    beam_kpc = BEAM_ARCSEC * (distance_mpc * 1.0e3) / 206265.0
    sigma_beam = ALPHA_BEAM * beam_kpc * v_obs / (r + beam_kpc)
    # asymmetric drift (global, type-dependent)
    frac = DRIFT_DWARF if galaxy_type == 'dwarf' else DRIFT_SPIRAL
    sigma_asym = frac * v_obs
    # outer-disk turbulence/warp proxy
    R_d = float(R_d_kpc) if R_d_kpc and R_d_kpc > 0 else 3.0
    shape = np.power(1.0 - np.exp(-r / R_d), P_TURB)
    sigma_turb = K_TURB * v_obs * shape
    return np.sqrt(sigma_base ** 2 + sigma_beam ** 2 + sigma_asym ** 2 + sigma_turb ** 2)


def v_baryon_effective(df, v_baryon_fallback: np.ndarray, ml_disk: float) -> np.ndarray:
    if hasattr(df, 'columns') and all(c in df.columns for c in ('vgas', 'vdisk', 'vbul')):
        vgas = np.asarray(df['vgas'].values, float)
        vdisk = np.asarray(df['vdisk'].values, float)
        vbul = np.asarray(df['vbul'].values, float)
        vdisk_scaled = vdisk * np.sqrt(max(ml_disk, 0.0))
        v2 = np.maximum(vgas, 0.0) ** 2 + np.maximum(vdisk_scaled, 0.0) ** 2 + np.maximum(vbul, 0.0) ** 2
        return np.sqrt(np.maximum(v2, 0.0))
    return np.asarray(v_baryon_fallback, float)


def evaluate_variant_time(master_table: Dict[str, Any], use_n_profile: bool, use_xi: bool, ml_disk: float) -> Dict[str, Any]:
    T_ref_years = collect_T_ref(master_table)
    results = []
    w_all = []
    for name, gal in master_table.items():
        df = gal.get("data")
        r_full = np.asarray(gal.get("r", []), dtype=float)
        v_obs_full = np.asarray(gal.get("v_obs", []), dtype=float)
        v_baryon_fb = np.asarray(gal.get("v_baryon", []), dtype=float)
        T_dyn_years = np.asarray(gal.get("T_dyn", []), dtype=float)
        R_d = float(gal.get("R_d", 3.0))
        if r_full.size < 2 or v_obs_full.size != r_full.size or v_baryon_fb.size != r_full.size or T_dyn_years.size != r_full.size:
            continue
        # Global mask: exclude points within beam scale
        beam_kpc = BEAM_ARCSEC * (float(gal.get('distance', 10.0)) * 1.0e3) / 206265.0
        mask = r_full >= beam_kpc
        r = r_full[mask]
        v_obs = v_obs_full[mask]
        v_baryon_fb = v_baryon_fb[mask]
        T_dyn_years = T_dyn_years[mask]
        if r.size < 2:
            continue
        # compute effective baryon using masked components
        if hasattr(df, 'columns') and all(c in df.columns for c in ('vgas','vdisk','vbul')):
            vgas = np.asarray(df['vgas'].values, float)[mask]
            vdisk = np.asarray(df['vdisk'].values, float)[mask]
            vbul = np.asarray(df['vbul'].values, float)[mask]
            vdisk_scaled = vdisk * np.sqrt(max(ml_disk, 0.0))
            v2 = np.maximum(vgas, 0.0) ** 2 + np.maximum(vdisk_scaled, 0.0) ** 2 + np.maximum(vbul, 0.0) ** 2
            v_baryon = np.sqrt(np.maximum(v2, 0.0))
        else:
            v_baryon = np.asarray(v_baryon_fb, float)
        v_max = float(np.max(v_obs)) if v_obs.size else 0.0
        galaxy_type = 'dwarf' if v_max < 80.0 else 'spiral'
        distance_mpc = float(gal.get('distance', 10.0))
        seff = sigma_components(df, r, v_obs, distance_mpc, galaxy_type, R_d, idx_mask=mask)
        w_t = w_small_lag_centered(T_dyn_years, T_ref_years)
        n_r = n_normalised(r) if use_n_profile else np.ones_like(r)
        zeta = zeta_vertical(r, R_d)
        xi = xi_global(gal) if use_xi else 1.0
        w_total = np.maximum(0.0, w_t) * np.maximum(0.0, n_r) * np.maximum(0.0, zeta) * float(xi)
        v_model = np.sqrt(w_total) * np.maximum(0.0, v_baryon)
        chi2 = float(np.sum(((v_obs - v_model) / seff) ** 2))
        w_all.append(w_total)
        results.append({
            "name": name,
            "chi2_reduced": chi2 / max(1, int(v_obs.size)),
        })
    chi_list = np.array([r["chi2_reduced"] for r in results]) if results else np.array([np.nan])
    w_cat = np.concatenate(w_all) if w_all else np.array([np.nan])
    summary = {
        "kernel": "time",
        "use_n_profile": bool(use_n_profile),
        "use_xi": bool(use_xi),
        "ml_disk": ml_disk,
        "c_lag": C_LAG,
        "A": A_N,
        "r0_kpc": R0_KPC,
        "p_shape": P_SHAPE,
        "R_weight_kpc": R_WEIGHT_KPC,
        "N_global": N_GLOBAL,
        "c_xi": C_XI,
        "gamma_xi": GAMMA_XI,
        "sigma0_kms": SIGMA0_KMS,
        "frac_floor": FRAC_FLOOR,
        "alpha_beam": ALPHA_BEAM,
        "drift_dwarf": DRIFT_DWARF,
        "drift_spiral": DRIFT_SPIRAL,
        "hz_ratio": HZ_RATIO,
        "N_galaxies": len(results),
        "chi2_reduced_median": float(np.nanmedian(chi_list)),
        "chi2_reduced_mean": float(np.nanmean(chi_list)),
        "w_p5": float(np.nanpercentile(w_cat, 5.0)),
        "w_p50": float(np.nanpercentile(w_cat, 50.0)),
        "w_p95": float(np.nanpercentile(w_cat, 95.0)),
        "top10": sorted(results, key=lambda r: r["chi2_reduced"], reverse=True)[:10],
    }
    return {"results": results, "summary": summary}


def evaluate_variant_accel(master_table: Dict[str, Any], use_n_profile: bool, use_xi: bool, ml_disk: float, g_ext_frac: float = 0.0) -> Dict[str, Any]:
    results = []
    w_all = []
    for name, gal in master_table.items():
        df = gal.get("data")
        r_full = np.asarray(gal.get("r", []), dtype=float)
        v_obs_full = np.asarray(gal.get("v_obs", []), dtype=float)
        v_baryon_fb = np.asarray(gal.get("v_baryon", []), dtype=float)
        T_dyn_years = np.asarray(gal.get("T_dyn", []), dtype=float)
        R_d = float(gal.get("R_d", 3.0))
        if r_full.size < 2 or v_obs_full.size != r_full.size or v_baryon_fb.size != r_full.size or T_dyn_years.size != r_full.size:
            continue
        beam_kpc = BEAM_ARCSEC * (float(gal.get('distance', 10.0)) * 1.0e3) / 206265.0
        mask = r_full >= beam_kpc
        r = r_full[mask]
        v_obs = v_obs_full[mask]
        v_baryon_fb = v_baryon_fb[mask]
        T_dyn_years = T_dyn_years[mask]
        if r.size < 2:
            continue
        if hasattr(df, 'columns') and all(c in df.columns for c in ('vgas','vdisk','vbul')):
            vgas = np.asarray(df['vgas'].values, float)[mask]
            vdisk = np.asarray(df['vdisk'].values, float)[mask]
            vbul = np.asarray(df['vbul'].values, float)[mask]
            vdisk_scaled = vdisk * np.sqrt(max(ml_disk, 0.0))
            v2 = np.maximum(vgas, 0.0) ** 2 + np.maximum(vdisk_scaled, 0.0) ** 2 + np.maximum(vbul, 0.0) ** 2
            v_baryon = np.sqrt(np.maximum(v2, 0.0))
        else:
            v_baryon = np.asarray(v_baryon_fb, float)
        v_max = float(np.max(v_obs)) if v_obs.size else 0.0
        galaxy_type = 'dwarf' if v_max < 80.0 else 'spiral'
        distance_mpc = float(gal.get('distance', 10.0))
        seff = sigma_components(df, r, v_obs, distance_mpc, galaxy_type, R_d, idx_mask=mask)
        w_g = w_accel_centered(v_baryon, r, g_ext_frac=g_ext_frac)
        n_r = n_normalised(r) if use_n_profile else np.ones_like(r)
        zeta = zeta_vertical(r, R_d)
        xi = xi_global(gal) if use_xi else 1.0
        w_total = np.maximum(0.0, w_g) * np.maximum(0.0, n_r) * np.maximum(0.0, zeta) * float(xi)
        v_model = np.sqrt(w_total) * np.maximum(0.0, v_baryon)
        chi2 = float(np.sum(((v_obs - v_model) / seff) ** 2))
        w_all.append(w_total)
        results.append({
            "name": name,
            "chi2_reduced": chi2 / max(1, int(v_obs.size)),
        })
    chi_list = np.array([r["chi2_reduced"] for r in results]) if results else np.array([np.nan])
    w_cat = np.concatenate(w_all) if w_all else np.array([np.nan])
    summary = {
        "kernel": "accel",
        "use_n_profile": bool(use_n_profile),
        "use_xi": bool(use_xi),
        "ml_disk": ml_disk,
        "c_lag": C_LAG,
        "A": A_N,
        "r0_kpc": R0_KPC,
        "p_shape": P_SHAPE,
        "R_weight_kpc": R_WEIGHT_KPC,
        "N_global": N_GLOBAL,
        "c_xi": C_XI,
        "gamma_xi": GAMMA_XI,
        "sigma0_kms": SIGMA0_KMS,
        "frac_floor": FRAC_FLOOR,
        "alpha_beam": ALPHA_BEAM,
        "drift_dwarf": DRIFT_DWARF,
        "drift_spiral": DRIFT_SPIRAL,
        "hz_ratio": HZ_RATIO,
        "a0_rs": A0_RS,
        "g_ext_frac": float(g_ext_frac),
        "N_galaxies": len(results),
        "chi2_reduced_median": float(np.nanmedian(chi_list)),
        "chi2_reduced_mean": float(np.nanmean(chi_list)),
        "w_p5": float(np.nanpercentile(w_cat, 5.0)),
        "w_p50": float(np.nanpercentile(w_cat, 50.0)),
        "w_p95": float(np.nanpercentile(w_cat, 95.0)),
        "top10": sorted(results, key=lambda r: r["chi2_reduced"], reverse=True)[:10],
    }
    return {"results": results, "summary": summary}


def evaluate_detailed(master_table: Dict[str, Any], kernel: str, use_n_profile: bool, use_xi: bool, ml_disk: float) -> Tuple[list, Dict[str, Any]]:
    """Compute per-galaxy details for a given configuration."""
    results = []
    w_all = []
    for name, gal in master_table.items():
        df = gal.get("data")
        r = np.asarray(gal.get("r", []), dtype=float)
        v_obs = np.asarray(gal.get("v_obs", []), dtype=float)
        v_baryon_fb = np.asarray(gal.get("v_baryon", []), dtype=float)
        T_dyn_years = np.asarray(gal.get("T_dyn", []), dtype=float)
        R_d = float(gal.get("R_d", 3.0))
        if r.size < 2 or v_obs.size != r.size or v_baryon_fb.size != r.size:
            continue
        v_baryon = v_baryon_effective(df, v_baryon_fb, ml_disk)
        v_max = float(np.max(v_obs)) if v_obs.size else 0.0
        galaxy_type = 'dwarf' if v_max < 80.0 else 'spiral'
        distance_mpc = float(gal.get('distance', 10.0))
        seff = sigma_components(df, r, v_obs, distance_mpc, galaxy_type, R_d)
        n_r = n_normalised(r) if use_n_profile else np.ones_like(r)
        zeta = zeta_vertical(r, R_d)
        xi = xi_global(gal) if use_xi else 1.0
        if kernel == 'time':
            T_ref_years = collect_T_ref(master_table)
            w_core = w_small_lag_centered(T_dyn_years, T_ref_years)
        else:
            w_core = w_accel_centered(v_baryon, r)
        w_total = np.maximum(0.0, w_core) * np.maximum(0.0, n_r) * np.maximum(0.0, zeta) * float(xi)
        v_model = np.sqrt(w_total) * np.maximum(0.0, v_baryon)
        resid = v_obs - v_model
        chi2 = float(np.sum(((resid) / seff) ** 2))
        chi2_red = chi2 / max(1, int(v_obs.size))
        w_all.append(w_total)
        results.append({
            "name": name,
            "chi2_reduced": chi2_red,
            "N": int(v_obs.size),
            "v_obs_mean": float(np.nanmean(v_obs)),
            "resid_rms": float(np.sqrt(np.nanmean(resid**2))),
        })
    chi_list = np.array([r["chi2_reduced"] for r in results]) if results else np.array([np.nan])
    w_cat = np.concatenate(w_all) if w_all else np.array([np.nan])
    summary = {
        "kernel": kernel,
        "use_n_profile": bool(use_n_profile),
        "use_xi": bool(use_xi),
        "ml_disk": ml_disk,
        "chi2_reduced_median": float(np.nanmedian(chi_list)),
        "chi2_reduced_mean": float(np.nanmean(chi_list)),
        "w_p50": float(np.nanpercentile(w_cat, 50.0)),
    }
    return results, summary


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    results_dir = root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    master_path = results_dir / "sparc_master.pkl"
    if not master_path.exists():
        print(f"ERROR: {master_path} not found. Build it with scripts/build_sparc_master_table.py")
        return
    with open(master_path, "rb") as f:
        master_table = pickle.load(f)

    # Evaluate across kernels, n-profile, xi, and ml_disk values
    bundle = {"time": {}, "accel": {}}
    rows_summary = []

    print("ILG pure solver summary (global variants, with σ floors, disk, beam/asym):")
    def pfx(tag, out):
        s = out["summary"]
        print(f"-- {tag} --")
        print(f"  kernel: {s['kernel']}  ml_disk: {s['ml_disk']}")
        print(f"  N_galaxies: {s['N_galaxies']}")
        print(f"  chi2_reduced_median: {s['chi2_reduced_median']}")
        print(f"  chi2_reduced_mean:   {s['chi2_reduced_mean']}")
        print(f"  w percentiles (p5, p50, p95): {s['w_p5']:.3f}, {s['w_p50']:.3f}, {s['w_p95']:.3f}")
        rows_summary.append([tag, s['kernel'], s.get('ml_disk',''), s['N_galaxies'], s['chi2_reduced_median'], s['chi2_reduced_mean']])

    for ml in ML_DISK_VALUES:
        out_t_with_n_no_xi = evaluate_variant_time(master_table, True, False, ml)
        out_t_with_n_xi    = evaluate_variant_time(master_table, True, True, ml)
        out_t_no_n_no_xi   = evaluate_variant_time(master_table, False, False, ml)
        out_t_no_n_xi      = evaluate_variant_time(master_table, False, True, ml)
        pfx(f"time: with n(r), no xi", out_t_with_n_no_xi)
        pfx(f"time: with n(r), xi",    out_t_with_n_xi)
        pfx(f"time: n(r)=1, no xi",    out_t_no_n_no_xi)
        pfx(f"time: n(r)=1, xi",       out_t_no_n_xi)
        bundle["time"][f"ml_{ml}"] = {
            "with_n_no_xi": out_t_with_n_no_xi,
            "with_n_xi": out_t_with_n_xi,
            "no_n_no_xi": out_t_no_n_no_xi,
            "no_n_xi": out_t_no_n_xi,
        }

        for gext in G_EXT_TESTS:
            out_g_with_n_no_xi = evaluate_variant_accel(master_table, True, False, ml, g_ext_frac=gext)
            out_g_with_n_xi    = evaluate_variant_accel(master_table, True, True, ml, g_ext_frac=gext)
            out_g_no_n_no_xi   = evaluate_variant_accel(master_table, False, False, ml, g_ext_frac=gext)
            out_g_no_n_xi      = evaluate_variant_accel(master_table, False, True, ml, g_ext_frac=gext)
            pfx(f"accel(gext={gext}): with n(r), no xi", out_g_with_n_no_xi)
            pfx(f"accel(gext={gext}): with n(r), xi",    out_g_with_n_xi)
            pfx(f"accel(gext={gext}): n(r)=1, no xi",    out_g_no_n_no_xi)
            pfx(f"accel(gext={gext}): n(r)=1, xi",       out_g_no_n_xi)
            bundle.setdefault("accel", {}).setdefault(f"ml_{ml}", {})[f"gext_{gext}"] = {
                "with_n_no_xi": out_g_with_n_no_xi,
                "with_n_xi": out_g_with_n_xi,
                "no_n_no_xi": out_g_no_n_no_xi,
                "no_n_xi": out_g_no_n_xi,
            }

    # Save pickle bundle
    out_path = results_dir / "ledger_final_pure_results.pkl"
    with open(out_path, "wb") as f:
        pickle.dump(bundle, f)
    print(f"Saved: {out_path}")

    # Write summary CSV
    summary_csv = results_dir / "ilg_variants_summary.csv"
    with open(summary_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["tag","kernel","ml_disk","N_gal","chi2_median","chi2_mean"])
        w.writerows(rows_summary)
    print(f"Saved: {summary_csv}")

    # Choose best median and best mean variants for detailed per-galaxy CSVs
    # Locked defaults: accel, with n, xi, ml=1.0
    best_median_results, best_median_summary = evaluate_detailed(master_table, kernel='accel', use_n_profile=True, use_xi=True, ml_disk=1.0)
    best_mean_results, best_mean_summary = evaluate_detailed(master_table, kernel='accel', use_n_profile=False, use_xi=True, ml_disk=1.0)

    def write_per_gal(fname: Path, results: list, summary: Dict[str, Any]):
        with open(fname, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(["kernel","use_n","use_xi","ml_disk","chi2_median","chi2_mean"]) 
            w.writerow([summary['kernel'], summary['use_n_profile'], summary['use_xi'], summary['ml_disk'], summary['chi2_reduced_median'], summary['chi2_reduced_mean']])
            w.writerow([])
            w.writerow(["name","N","chi2_reduced","v_obs_mean","resid_rms"]) 
            for row in sorted(results, key=lambda r: r['chi2_reduced'], reverse=True):
                w.writerow([row['name'], row['N'], row['chi2_reduced'], row['v_obs_mean'], row['resid_rms']])
        print(f"Saved: {fname}")

    write_per_gal(results_dir / "ilg_best_median_per_galaxy.csv", best_median_results, best_median_summary)
    write_per_gal(results_dir / "ilg_best_mean_per_galaxy.csv", best_mean_results, best_mean_summary)

    # Optional ablation mode (env var or flag-free trigger by default last step)
    # Minimal internal ablation: stepwise enabling of n, xi, drift, mask, zcap, turb
    steps = [
        ("base",              dict(n=False, xi=False, ml=1.0, turb=False, drift=False, mask=False, zcap=False)),
        ("+ n(r)",            dict(n=True,  xi=False, ml=1.0, turb=False, drift=False, mask=False, zcap=False)),
        ("+ xi",              dict(n=True,  xi=True,  ml=1.0, turb=False, drift=False, mask=False, zcap=False)),
        ("+ drift",           dict(n=True,  xi=True,  ml=1.0, turb=False, drift=True,  mask=False, zcap=False)),
        ("+ mask",            dict(n=True,  xi=True,  ml=1.0, turb=False, drift=True,  mask=True,  zcap=False)),
        ("+ zcap",            dict(n=True,  xi=True,  ml=1.0, turb=False, drift=True,  mask=True,  zcap=True)),
        ("+ turb (k=0.07,p=1.3)", dict(n=True,  xi=True,  ml=1.0, turb=True,  drift=True,  mask=True,  zcap=True)),
    ]
    # Save and restore globals used in toggles
    K_old, P_old = K_TURB, P_TURB
    DD_old, DS_old = DRIFT_DWARF, DRIFT_SPIRAL
    BA_old = BEAM_ARCSEC
    ZMIN_old, ZMAX_old = ZCAP_MIN, ZCAP_MAX

    ablation_rows = []
    for tag, cfg in steps:
        # apply toggles
        # turbulence
        globals()["K_TURB"], globals()["P_TURB"] = (0.0, 1.0) if not cfg["turb"] else (0.07, 1.3)
        # drift
        globals()["DRIFT_DWARF"], globals()["DRIFT_SPIRAL"] = ((0.0, 0.0) if not cfg["drift"] else (0.10, 0.05))
        # mask
        globals()["BEAM_ARCSEC"] = (0.0 if not cfg["mask"] else 15.0)
        # zcap
        globals()["ZCAP_MIN"], globals()["ZCAP_MAX"] = ((0.5, 1.5) if not cfg["zcap"] else (0.8, 1.2))
        # evaluate best kernel (accel)
        out = evaluate_variant_accel(master_table, use_n_profile=cfg["n"], use_xi=cfg["xi"], ml_disk=cfg["ml"], g_ext_frac=0.0)
        s = out["summary"]
        ablation_rows.append([tag, s["chi2_reduced_median"], s["chi2_reduced_mean"]])
        print(f"Ablation {tag:20s} -> median {s['chi2_reduced_median']:.3f}, mean {s['chi2_reduced_mean']:.3f}")

    # restore
    globals()["K_TURB"], globals()["P_TURB"] = K_old, P_old
    globals()["DRIFT_DWARF"], globals()["DRIFT_SPIRAL"] = DD_old, DS_old
    globals()["BEAM_ARCSEC"] = BA_old
    globals()["ZCAP_MIN"], globals()["ZCAP_MAX"] = ZMIN_old, ZMAX_old

    # write ablation csv
    abl_csv = results_dir / "ilg_ablation.csv"
    with open(abl_csv, 'w', newline='') as f:
        import csv
        w = csv.writer(f)
        w.writerow(["step","chi2_median","chi2_mean"])
        w.writerows(ablation_rows)
    print(f"Saved: {abl_csv}")

    # Figures: (1) v_obs_mean vs resid_rms (dwarf vs spiral) for locked best config
    #          (2) chi2 hist dwarf vs spiral
    #          (3) Ablation bar of chi2_median
    # Best config details
    best_res = best_median_results
    names = [r['name'] for r in best_res]
    vmeans = np.array([r['v_obs_mean'] for r in best_res])
    resid = np.array([r['resid_rms'] for r in best_res])
    dwarfs = vmeans < 80.0
    # (1) scatter
    plt.figure(figsize=(6,4))
    plt.scatter(vmeans[dwarfs], resid[dwarfs], s=18, c='#1f77b4', label='dwarfs')
    plt.scatter(vmeans[~dwarfs], resid[~dwarfs], s=18, c='#ff7f0e', label='spirals')
    plt.xlabel('Mean observed velocity [km/s]')
    plt.ylabel('Residual RMS [km/s]')
    plt.grid(alpha=0.3)
    plt.legend()
    fig1 = results_dir / 'ilg_resid_vs_vmean.png'
    plt.tight_layout(); plt.savefig(fig1, dpi=160)
    print(f"Saved: {fig1}")
    plt.close()
    # (2) chi2 hist dwarf vs spiral
    chis = np.array([r['chi2_reduced'] for r in best_res])
    plt.figure(figsize=(6,4))
    plt.hist(chis[dwarfs], bins=20, alpha=0.6, label='dwarfs')
    plt.hist(chis[~dwarfs], bins=20, alpha=0.6, label='spirals')
    plt.xlabel('Reduced chi-squared')
    plt.ylabel('Count')
    plt.grid(alpha=0.3)
    plt.legend()
    fig2 = results_dir / 'ilg_chi2_dwarf_spiral.png'
    plt.tight_layout(); plt.savefig(fig2, dpi=160)
    print(f"Saved: {fig2}")
    plt.close()
    # (3) Ablation bar
    steps_lbl = [row[0] for row in ablation_rows]
    chi_med = [row[1] for row in ablation_rows]
    plt.figure(figsize=(8,4))
    plt.bar(range(len(steps_lbl)), chi_med, color='#2ca02c')
    plt.xticks(range(len(steps_lbl)), steps_lbl, rotation=30, ha='right')
    plt.ylabel('Median chi-squared/N')
    plt.grid(axis='y', alpha=0.3)
    fig3 = results_dir / 'ilg_ablation_bar.png'
    plt.tight_layout(); plt.savefig(fig3, dpi=160)
    print(f"Saved: {fig3}")
    plt.close()


if __name__ == "__main__":
    main()
