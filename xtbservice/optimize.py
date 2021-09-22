
import io
from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
from functools import lru_cache

from ase import Atoms
from ase.optimize.lbfgs import LBFGS
from xtb.ase.calculator import XTB

from .cache import opt_cache
from .models import OptimizationResult
from .utils import hash_atoms


def opt_hash(atoms, method):
    return hash(str(hash_atoms(atoms)) + method)

def run_xtb_opt(
    atoms: Atoms, method: str="GFNFF", constrain_metals: bool = False, fmax: float = 0.00001
) -> OptimizationResult:
    this_hash = opt_hash(atoms, method)
    try:
        result = opt_cache.get(this_hash)
    except KeyError:
        pass
    if result is None:
        mol = deepcopy(atoms)
        mol.pbc = False
        

        mol.calc = XTB(method=method)
        opt = LBFGS(mol, logfile=None)
        opt.run(fmax=fmax)
        forces = mol.get_forces()
        energy = mol.get_potential_energy()
        mol.calc = None
        opt_cache.set(this_hash, result,  expire=None)
        result =  OptimizationResult(atoms=mol, forces=forces, energy=energy)
    return result
