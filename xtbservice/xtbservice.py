# -*- coding: utf-8 -*-
"""
xtb-service.py
webservice providing xtb calculations
"""
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi_versioning import VersionedFastAPI, version
from starlette.middleware import Middleware

from . import __version__
from .conformers import conformers_from_molfile, conformers_from_smiles
from .ir import ir_from_molfile, ir_from_smiles
from .models import ConformerLibrary, ConformerRequest, IRRequest, IRResult

ALLOWED_HOSTS = ["*"]


app = FastAPI(
    title="XTB webservice",
    description="Offers xtb calculation tools. Allowed methods are `GFNFF`, `GFN2xTB`, `GFN1xTB`",
    version=__version__,
    contact={
        "name": "Cheminfo",
        "email": "admin@cheminfo.org",
    },
    license_info={"name": "MIT"},
)


@app.get("/app_version")
@version(1)
def read_version():
    return {"app_version": __version__}


@app.post("/ir", response_model=IRResult)
@version(1)
def post_get_ir_spectrum(irrequest: IRRequest):
    if irrequest.smiles:
        ir = ir_from_smiles(irrequest.smiles, irrequest.method)
    elif irrequest.molFile:
        ir = ir_from_molfile(irrequest.molFile, irrequest.method)
    else:
        raise HTTPException(
            status_code=422, detail="You need to provide either `molFile` or `smiles`"
        )
    return ir


@app.post("/conformers", response_model=ConformerLibrary)
@version(1)
def post_conformers(conformerrequest: ConformerRequest):
    if conformerrequest.smiles:
        conformers = conformers_from_smiles(
            conformerrequest.smiles,
            conformerrequest.forceField,
            conformerrequest.rmsdThreshold,
            conformerrequest.maxConformers,
        )
    elif conformerrequest.molFile:
        conformers = conformers_from_molfile(
            conformerrequest.molFile,
            conformerrequest.forceField,
            conformerrequest.rmsdThreshold,
            conformerrequest.maxConformers,
        )
    else:
        raise HTTPException(
            status_code=422, detail="You need to provide either `molFile` or `smiles`"
        )
    return conformers


@app.get("/ir", response_model=IRResult)
@version(1)
def get_ir_spectrum(smiles: str, method: str = "GFNFF"):
    ir = ir_from_smiles(smiles, method)
    return ir


app = VersionedFastAPI(
    app,
    version_format="{major}",
    prefix_format="/v{major}",
    middleware=[
        Middleware(
            CORSMiddleware,
            allow_origins=ALLOWED_HOSTS,
            allow_credentials=True,
            allow_methods=["*"],
            allow_headers=["*"],
        )
    ],
)
