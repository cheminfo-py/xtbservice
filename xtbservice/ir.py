from .models import IRResult
from copy import deepcopy
from ase.vibrations import Infrared
from xtb.ase.calculator import XTB
from ase import Atoms
from .utils import smiles2ase, hash_atoms
from .optimize import run_xtb_opt
from functools import lru_cache
from .cache import ir_cache
import numpy as np


def ir_hash(atoms, method):
    return hash(str(hash_atoms(atoms)) + method)


def run_xtb_ir(atoms: Atoms, method: str = "GFNFF") -> IRResult:
    # mol = deepcopy(atoms)
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
        spectrum = ir.get_spectrum()
        zpe = ir.get_zero_point_energy()
        ir.clean()
        result = IRResult(
            wavenumbers=list(spectrum[0]),
            intensities=list(spectrum[1]),
            zeroPointEnergy=zpe,
            displacements=get_displacement_xyz_dict(ir),
            mostRelevantModesOfAtoms=get_max_displacements(ir),
        )
        ir_cache.set(this_hash, result)
    return result


@lru_cache()
def ir_from_smiles(smiles, method):
    atoms = smiles2ase(smiles)
    opt_result = run_xtb_opt(atoms, method=method)
    ir = run_xtb_ir(opt_result.atoms, method=method)
    return ir



def get_max_displacements(ir): 
    mode_abs_displacements = []
    
    for n in range(3 * len(ir.indices)):
        mode_abs_displacements.append(np.linalg.norm(ir.get_mode(n), axis=1)) 
    
    mode_abs_displacements = np.stack(mode_abs_displacements)
    return dict(zip(ir.indices, [list(a) for a in mode_abs_displacements.argsort(axis=0)]))


def get_displacement_xyz_dict(ir):
    symbols = ir.atoms.get_chemical_symbols()
    freq = ir.get_frequencies()
    modes = {}

    for n in range(3 * len(ir.indices)):
        xyz_file = []
        xyz_file.append("%6d\n" % len(ir.atoms))

        if freq[n].imag != 0:
            c = "i"
            freq[n] = freq[n].imag

        else:
            freq[n] = freq[n].real
            c = " "

        xyz_file.append("Mode #%d, f = %.1f%s cm^-1" % (n, float(freq[n].real), c))

        if ir.ir:
            xyz_file.append(", I = %.4f (D/Å)^2 amu^-1.\n" % ir.intensities[n])
        else:
            xyz_file.append(".\n")

        # dict_label = xyz_file[-1] + xyz_file[-2]
        # dict_label = dict_label.strip('\n')

        mode = ir.get_mode(n)
        for i, pos in enumerate(ir.atoms.positions):
            xyz_file.append(
                "%2s %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n"
                % (
                    symbols[i],
                    pos[0],
                    pos[1],
                    pos[2],
                    mode[i, 0],
                    mode[i, 1],
                    mode[i, 2],
                )
            )

        xyz_file_string = "".join(xyz_file)

        modes[n] = xyz_file_string

    return modes
