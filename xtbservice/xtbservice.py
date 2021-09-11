"""
xtb-service.py
webservice providing xtb calculations
"""
from . import __version__
from fastapi import FastAPI
from .ir import run_xtb_ir
from .utils import smiles2ase
from .optimize import run_xtb_optimize
from .models import IRRequest

app = FastAPI()


@app.get("/version")
def read_version():
    return {"version": __version__}

@app.get("/ir")
def get_ir_spectrum(IRRequest):
    atoms = smiles2ase(IRRequest.smiles)
    atoms = run_xtb_optimize(atoms, method=IRRequest.method)
    ir = run_xtb_ir(atoms, method=IRRequest.method)
    return ir
