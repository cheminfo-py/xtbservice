from dataclasses import dataclass
import numpy as np 
from pydantic import BaseModel
from typing import Optional

@dataclass
class OptimizationResult:
    atoms: Atoms
    forces: np.ndarray
    energy: float


@dataclass
class IRResult:
    wavenumber: np.array
    intensity: np.array
    zero_point_energy: float

class IRRequest(BaseModel):
    smiles: str
    method: Optional[str] = 'GFNFF'