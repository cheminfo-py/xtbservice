import io
import shutil
from contextlib import redirect_stdout
from functools import lru_cache
from rdkit import Chem
from typing import Union, List
import numpy as np
from ase import Atoms
from ase.vibrations import Infrared
from xtb.ase.calculator import XTB

from .cache import ir_cache, ir_from_molfile_cache, ir_from_smiles_cache
from .models import IRResult
from .optimize import run_xtb_opt
from .utils import get_hash, hash_atoms, molfile2ase, smiles2ase


def ir_hash(atoms, method):
    return hash(str(hash_atoms(atoms)) + method)


def run_xtb_ir(
    atoms: Atoms, method: str = "GFNFF", mol: Union[None, Chem.Mol] = None
) -> IRResult:
    # mol = deepcopy(atoms)
    this_hash = ir_hash(atoms, method)

    result = ir_cache.get(this_hash)

    if result is None:
        atoms.pbc = False
        atoms.calc = XTB(method=method)

        ir = Infrared(atoms, name=str(this_hash))
        ir.run()
        spectrum = ir.get_spectrum(start=500, end=4000)
        zpe = ir.get_zero_point_energy()
        most_relevant_mode_for_bond = None
        bond_displacements = None
        if mol is not None:
            bond_displacements = compile_all_bond_displacements(mol, atoms, ir)
            most_relevant_mode_for_bond_ = np.argmax(bond_displacements, axis=0)
            most_relevant_mode_for_bond = []
            bonds = get_bonds_from_mol(mol)
            for i, mode in enumerate(most_relevant_mode_for_bond_):
                most_relevant_mode_for_bond.append(
                    {"startAtom": bonds[i][0], "endAtom": bonds[i][1], "mode": int(mode), "displacement": bond_displacements[mode][i]}
                )
         
        mode_info, has_imaginary = compile_modes_info(ir, bond_displacements, bonds)
        result = IRResult(
            wavenumbers=list(spectrum[0]),
            intensities=list(spectrum[1]),
            zeroPointEnergy=zpe,
            modes=mode_info,
            hasImaginaryFrequency=has_imaginary,
            mostRelevantModesOfAtoms=get_max_displacements(ir),
            mostRelevantModesOfBonds=most_relevant_mode_for_bond
        )
        ir_cache.set(this_hash, result)

        shutil.rmtree(ir.cache.directory)
        ir.clean()
    return result


def ir_from_smiles(smiles, method):
    myhash = str(get_hash(smiles + method))

    result = ir_from_smiles_cache.get(myhash)

    if result is None:
        atoms, mol = smiles2ase(smiles)
        opt_result = run_xtb_opt(atoms, method=method)
        result = run_xtb_ir(opt_result.atoms, method=method, mol=mol)
        ir_from_smiles_cache.set(myhash, result, expire=None)
    return result


def ir_from_molfile(molfile, method):
    myhash = str(get_hash(molfile + method))

    result = ir_from_molfile_cache.get(myhash)

    if result is None:
        atoms, mol = molfile2ase(molfile)
        opt_result = run_xtb_opt(atoms, method=method)
        result = run_xtb_ir(opt_result.atoms, method=method, mol=mol)
        ir_from_molfile_cache.set(myhash, result, expire=None)
    return result


def compile_all_bond_displacements(mol, atoms, ir):
    bond_displacements = []
    for mode_number in range(3 * len(ir.indices)):
        bond_displacements.append(
            get_bond_displacements(mol, atoms, ir.get_mode(mode_number))
        )

    return np.vstack(bond_displacements)


def clean_frequency(frequencies, n):
    if frequencies[n].imag != 0:
        c = "i"
        freq = frequencies[n].imag

    else:
        freq = frequencies[n].real
        c = " "
    return freq, c


def compile_modes_info(ir, bond_displacements=None, bonds=None):
    frequencies = ir.get_frequencies()
    symbols = ir.atoms.get_chemical_symbols()
    modes = []
    has_imaginary = False
    for n in range(3 * len(ir.indices)):
        f, c = clean_frequency(frequencies, n)
        has_imaginary = True if c == "i" else False
        mostContributingBonds = None
        if bond_displacements is not None: 
            mostContributingBonds = select_most_contributing_bonds(bond_displacements[n,:])
            mostContributingBonds = [bonds[i] for i in mostContributingBonds]
        modes.append(
            {
                "number": n,
                "displacements": get_displacement_xyz_for_mode(
                    ir, frequencies, symbols, n
                ),
                "intensity": float(ir.intensities[n]),
                "frequency": float(f),
                "imaginary": True if c == "i" else False,
                "mostDisplacedAtoms": [
                    int(i) for i in np.argsort(np.linalg.norm(ir.get_mode(n), axis=1))
                ][::-1],
                "mostContributingAtoms": [
                    int(i) for i in select_most_contributing_atoms(ir, n)
                ],
                "mostContributingBonds": mostContributingBonds
            }
        )

    return modes, has_imaginary


def get_max_displacements(ir):
    mode_abs_displacements = []

    for n in range(3 * len(ir.indices)):
        mode_abs_displacements.append(np.linalg.norm(ir.get_mode(n), axis=1))

    mode_abs_displacements = np.stack(mode_abs_displacements)
    return dict(
        zip(ir.indices, [list(a) for a in mode_abs_displacements.argsort(axis=0)][::-1])
    )


def get_displacement_xyz_for_mode(ir, frequencies, symbols, n):
    xyz_file = []
    xyz_file.append("%6d\n" % len(ir.atoms))

    f, c = clean_frequency(frequencies, n)

    xyz_file.append("Mode #%d, f = %.1f%s cm^-1" % (n, float(f.real), c))

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
            % (symbols[i], pos[0], pos[1], pos[2], mode[i, 0], mode[i, 1], mode[i, 2],)
        )

    xyz_file_string = "".join(xyz_file)
    return xyz_file_string


def get_displacement_xyz_dict(ir):
    symbols = ir.atoms.get_chemical_symbols()
    frequencies = ir.get_frequencies()
    modes = {}

    for n in range(3 * len(ir.indices)):
        modes[n] = get_displacement_xyz_for_mode(ir, frequencies, symbols, n)

    return modes


def select_most_contributing_atoms(ir, mode, threshold: float = 0.4):
    relative_contribution = (
        np.linalg.norm(ir.get_mode(mode), axis=1)
        / np.linalg.norm(ir.get_mode(mode), axis=1).sum()
    )

    return np.where(
        relative_contribution
        > threshold * np.max(np.abs(np.diff(relative_contribution)))
    )[0]



def select_most_contributing_bonds(displacements, threshold: float = 0.4):
    relative_contribution = displacements / displacements.sum()

    return np.where(
        relative_contribution
        > threshold * np.max(np.abs(np.diff(relative_contribution)))
    )[0]


@lru_cache()
def get_bonds_from_mol(mol) -> List[tuple]:
    all_bonds = []
    for i in range(mol.GetNumBonds()):
        bond = mol.GetBondWithIdx(i)
        all_bonds.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

    return all_bonds


def get_bond_vector(positions, bond):
    return positions[bond[1]] - positions[bond[0]]


def get_displaced_positions(positions, mode):
    return positions + mode


def get_bond_displacements(mol, atoms, mode):
    bonds = get_bonds_from_mol(mol)
    positions = atoms.positions
    displaced_positions = get_displaced_positions(positions, mode)
    changes = []

    for bond in bonds:
        changes.append(
            np.linalg.norm(
                get_bond_vector(positions, bond)
                - get_bond_vector(displaced_positions, bond)
            )
        )

    return changes

