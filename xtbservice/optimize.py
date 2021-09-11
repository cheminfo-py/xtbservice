
from xtb.ase.calculator import XTB
from ase.optimize.lbfgs import LBFGS
from copy import deepcopy
from ase import Atoms

from .models import OptimizationResult


def run_xtb_opt(
    atoms: Atoms, method: str="GFNFF", constrain_metals: bool = False, fmax: float = 0.0001
) -> OptimizationResult:
    mol = deepcopy(atoms)
    mol.pbc = False
    mol.calc = XTB(method=method)
    opt = LBFGS(mol, logfile=None)
    opt.run(fmax=fmax)
    forces = mol.get_forces()
    energy = mol.get_potential_energy()
    return OptimizationResult(atoms=mol, forces=forces, energy=energy)


