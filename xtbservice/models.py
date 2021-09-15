from dataclasses import dataclass
import numpy as np
from pydantic import BaseModel, validator
from typing import Optional, List, Dict
from ase import Atoms

ALLOWED_METHODS = ('GFNFF', 'GFN2xTB', 'GFN1xTB')

@dataclass
class OptimizationResult:
    atoms: Atoms
    forces: np.ndarray
    energy: float


class IRResult(BaseModel):
    wavenumbers: List[float]
    intensities: List[float]
    zeroPointEnergy: float
    modes: Optional[List[dict]]
    mostRelevantModesOfAtoms: Optional[Dict[int, List[int]]]


class IRRequest(BaseModel):
    smiles: str
    method: Optional[str] = "GFNFF"

    @validator('method')
    def method_match(cls, v):
        if not v in ALLOWED_METHODS:
            raise ValueError(f'method must be in {ALLOWED_METHODS}')
        return v