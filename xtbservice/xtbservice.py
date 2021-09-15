"""
xtb-service.py
webservice providing xtb calculations
"""
from . import __version__
from fastapi import FastAPI
from .ir import ir_from_smiles
from .models import IRRequest, IRResult

app = FastAPI()


@app.get("/version")
def read_version():
    return {"version": __version__}


@app.post("/ir", response_model=IRResult)
def post_get_ir_spectrum(irrequest: IRRequest):
    ir = ir_from_smiles(irrequest.smiles, irrequest.method)
    return ir



@app.get("/ir", response_model=IRResult)
def get_ir_spectrum(smiles: str, method: str = 'GFNFF'):
    ir = ir_from_smiles(smiles, method)
    return ir
