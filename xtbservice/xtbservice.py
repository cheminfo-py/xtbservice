"""
xtb-service.py
webservice providing xtb calculations
"""
from . import __version__
from fastapi import FastAPI
from .ir import ir_from_smiles
from .models import IRRequest, IRResult
from fastapi.middleware.cors import CORSMiddleware

ALLOWED_HOSTS = ["*"]


app = FastAPI(
    title="XTB webservice",
    description="Offers xtb calculation tools. Allowed methods are `GFNFF`, `GFN2xTB`, `GFN1xTB`",
    version=__version__,
    contact={
        "name": "Cheminfo",
        "email": "admin@cheminfo.org",
    },
    license_info={
        "name": "MIT"}
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_HOSTS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


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
