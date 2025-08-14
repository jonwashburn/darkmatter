import importlib
import sys

# Avoid heavy/stochastic ML frameworks
FORBIDDEN = {"torch", "tensorflow", "jax"}


def test_no_forbidden_imports():
    modname = "scripts.ledger_final_combined"
    importlib.import_module(modname)
    for name, module in list(sys.modules.items()):
        if module is None:
            continue
        base = name.split(".")[0]
        assert base not in FORBIDDEN, f"Forbidden module imported: {base} in {modname}"
