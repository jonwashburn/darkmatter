# darkmatter: ILG rotation curves (Recognition Science)

This repository contains code and artifacts supporting the Information‑Limited Gravity (ILG) rotation‑curve analysis.

## Quick start

```bash
# build
docker build -t ilg-validation .
# run (uses bundled processed SPARC master table)
docker run --rm -v "$PWD":/app ilg-validation python scripts/ledger_final_combined.py | cat
```

## Environment (without Docker)

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/ledger_final_combined.py | cat
```

## Contents
- scripts/ledger_final_combined.py — main non‑relativistic solver
- scripts/build_sparc_master_table.py — create `results/sparc_master.pkl` from raw SPARC
- relativistic_rs_gravity.py — minimal prospective relativistic helpers
- results/ — processed data and outputs (includes `sparc_master.pkl`)
- figures/ — example plots
- RecognitionScience/, lean_proofs/ — derivation and Lean artifacts

## Tests

```bash
pytest -q | cat
```

## Notes
- Raw SPARC files are not included; the pipeline expects them if you want to rebuild `sparc_master.pkl`.
- The repository avoids per‑galaxy tuning; constants are fixed globally.
