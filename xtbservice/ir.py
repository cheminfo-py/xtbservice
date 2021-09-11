from .models import IRResult
from copy import deepcopy
from ase.vibrations import Infrared
from xtb.ase.calculator import XTB
from ase import Atoms
from .utils import smiles2ase
from .optimize import run_xtb_opt
from functools import lru_cache

def run_xtb_ir(atoms: Atoms, method: str="GFNFF") -> IRResult:
    #mol = deepcopy(atoms)
    atoms.pbc = False
    atoms.calc = XTB(method=method)
    ir = Infrared(atoms)
    ir.run()
    spectrum =  ir.get_spectrum()
    zpe = ir.get_zero_point_energy()
    return IRResult(
        wavenumber = list(spectrum[0]),
        intensities = list(spectrum[1]),
        zero_point_energy=zpe
    )

@lru_cache()
def ir_from_smiles(smiles, method):
    atoms = smiles2ase(smiles)
    atoms = run_xtb_opt(atoms, method=method)
    ir = run_xtb_ir(atoms.atoms, method=method)
    return ir

