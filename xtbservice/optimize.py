# -*- coding: utf-8 -*-


from copy import deepcopy

from ase import Atoms
from ase.optimize.lbfgs import LBFGS
from fastapi.logger import logger
from xtb.ase.calculator import XTB

from .cache import opt_cache
from .models import OptimizationResult
from .utils import hash_atoms, hash_object


def opt_hash(atoms, method):
    return hash_object(str(hash_atoms(atoms)) + method)


def run_xtb_opt(
    atoms: Atoms,
    method: str = "GFNFF",
    constrain_metals: bool = False,
    fmax: float = 0.000005,
    maxiter: int = 100,
) -> OptimizationResult:
    this_hash = opt_hash(atoms, method)
    logger.debug(f"Running optimization with hash {this_hash}")
    try:
        result = opt_cache.get(this_hash)
    except KeyError:
        pass
    if result is None:
        logger.debug(f"Optimization not found in cache, running")
        mol = deepcopy(atoms)
        mol.pbc = False

        mol.calc = XTB(method=method)
        opt = LBFGS(mol, logfile=None)
        opt.run(fmax=fmax, steps=maxiter)
        forces = mol.get_forces()
        energy = mol.get_potential_energy()
        mol.calc = None
        opt_cache.set(this_hash, result, expire=None)
        result = OptimizationResult(atoms=mol, forces=forces, energy=energy)
    return result
