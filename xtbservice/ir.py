from .models import IRResult
from copy import deepcopy
from ase.vibrations import Infrared
from xtb.ase.calculator import XTB
from ase import Atoms
from .utils import smiles2ase, hash_atoms
from .optimize import run_xtb_opt
from functools import lru_cache
from .cache import ir_cache

def ir_hash(atoms, method):
    return hash(str(hash_atoms(atoms)) + method)

def run_xtb_ir(atoms: Atoms, method: str="GFNFF") -> IRResult:
    #mol = deepcopy(atoms)
    this_hash = ir_hash(atoms, method)
    try:
        result = ir_cache.get(this_hash)
    except KeyError:
        pass
    if result is None:
        atoms.pbc = False
        atoms.calc = XTB(method=method)
        ir = Infrared(atoms, name=str(this_hash))
        ir.run()
        spectrum =  ir.get_spectrum()
        zpe = ir.get_zero_point_energy()
        ir.clean()
        result =  IRResult(
            wavenumbers = list(spectrum[0]),
            intensities = list(spectrum[1]),
            zero_point_energy=zpe
        )
        #ir_cache.set(this_hash, result)
    return result

@lru_cache()
def ir_from_smiles(smiles, method):
    atoms = smiles2ase(smiles)
    opt_result = run_xtb_opt(atoms, method=method)
    ir = run_xtb_ir(opt_result.atoms, method=method)
    return ir

