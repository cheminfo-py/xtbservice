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
from .errors import TooLargeError
from .ir import ir_from_molfile, ir_from_smiles
from .models import ConformerLibrary, ConformerRequest, IRRequest, IRResult
from .settings import MAX_ATOMS_FF, MAX_ATOMS_XTB

ALLOWED_HOSTS = ["*"]


app = FastAPI(
    title="XTB webservice",
    description="Offers xtb calculation tools. Allowed methods are `GFNFF`, `GFN2xTB`, `GFN1xTB`",
    version=__version__,
    contact={"name": "Cheminfo", "email": "admin@cheminfo.org",},
    license_info={"name": "MIT"},
)


def max_atoms_error():
    return HTTPException(
        status_code=422,
        detail=f"This services only accepts structures with less than {MAX_ATOMS_FF} atoms for force-field calculations and {MAX_ATOMS_XTB} for xtb calculations.",
    )


@app.get("/app_version")
@version(1)
def read_version():
    return {"app_version": __version__}


@app.post("/ir", response_model=IRResult)
@version(1)
def post_get_ir_spectrum(irrequest: IRRequest):
    try:
        if irrequest.smiles:
            ir = ir_from_smiles(irrequest.smiles, irrequest.method)
        elif irrequest.molFile:
            ir = ir_from_molfile(irrequest.molFile, irrequest.method)
        else:
            raise HTTPException(
                status_code=422,
                detail="You need to provide either `molFile` or `smiles`",
            )
    except TooLargeError:
        raise max_atoms_error()
    except TimeoutError:
        raise HTTPException(status_code=500, detail="Calculation timed out.")
    return ir


@app.post("/conformers", response_model=ConformerLibrary)
@version(1)
def post_conformers(conformerrequest: ConformerRequest):
    try:
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
                status_code=422,
                detail="You need to provide either `molFile` or `smiles`",
            )
    except TooLargeError:
        raise max_atoms_error()
    except TimeoutError:
        raise HTTPException(status_code=500, detail="Calculation timed out.")
    return conformers


@app.get("/ir", response_model=IRResult)
@version(1)
def get_ir_spectrum(smiles: str, method: str = "GFNFF"):
    try:
        ir = ir_from_smiles(smiles, method)
    except TooLargeError:
        raise max_atoms_error()
    except TimeoutError:
        raise HTTPException(status_code=500, detail="Calculation timed out.")
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
