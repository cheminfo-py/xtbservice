from .models import IRResult
from copy import deepcopy
from ase.vibrations import Infrared
from xtb.ase.calculator import XTB
from ase import Atoms

def run_xtb_ir(atoms: Atoms, method: str="GFNFF") -> IRResult:
    mol = deepcopy(atoms)
    mol.pbc = False
    mol.calc = XTB(method=method)
    ir = Infrared(atoms)
    ir.run()
    spectrum =  ir.get_spectrum()
    zpe = spectrum.get_zpe()
    return IRResult(
        wavenumber = spectrum[0],
        intensities = spectrum[1],
        zero_point_energy=zpe
    )

