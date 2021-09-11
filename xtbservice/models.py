from dataclasses import dataclass
import numpy as np
from pydantic import BaseModel
from typing import Optional, List
from ase import Atoms


@dataclass
class OptimizationResult:
    atoms: Atoms
    forces: np.ndarray
    energy: float


class IRResult(BaseModel):
    wavenumber: List[float]
    intensities: List[float]
    zero_point_energy: float


class IRRequest(BaseModel):
    smiles: str
    method: Optional[str] = "GFNFF"

