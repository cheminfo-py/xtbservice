# -*- coding: utf-8 -*-
import hashlib

import numpy as np
from ase import Atoms
from rdkit import Chem

from .cache import conformer_cache
from .conformers import embed_conformer
from .errors import TooLargeError
from .settings import MAX_ATOMS


def rdkit2ase(mol):
    pos = mol.GetConformer().GetPositions()
    natoms = mol.GetNumAtoms()
    species = [mol.GetAtomWithIdx(j).GetSymbol() for j in range(natoms)]
    atoms = Atoms(species, positions=pos)
    atoms.pbc = False

    return atoms


def molfile2ase(molfile: str) -> Atoms:
    try:
        result = conformer_cache.get(molfile)
    except KeyError:
        pass

    if result is None:
        mol = Chem.MolFromMolBlock(molfile, sanitize=True, removeHs=False)
        mol.UpdatePropertyCache(strict=False)
        natoms = mol.GetNumAtoms()
        if natoms > MAX_ATOMS:
            raise TooLargeError(
                f"Molecule can have maximal {MAX_ATOMS} atoms for this service"
            )
        mol = embed_conformer(mol)
        result = rdkit2ase(mol), mol
        conformer_cache.set(molfile, result, expire=None)
    return result


def smiles2ase(smiles: str) -> Atoms:
    try:
        result = conformer_cache.get(smiles)
    except KeyError:
        pass

    if result is None:
        mol = Chem.MolFromSmiles(smiles)
        natoms = mol.GetNumAtoms()
        if natoms > MAX_ATOMS:
            raise TooLargeError(
                f"Molecule can have maximal {MAX_ATOMS} atoms for this service"
            )
        refmol = Chem.AddHs(Chem.Mol(mol))
        refmol = embed_conformer(refmol)
        result = rdkit2ase(refmol), refmol
        conformer_cache.set(smiles, result, expire=None)
    return result


def hash_atoms(atoms: Atoms) -> int:
    symbols = str(atoms.symbols)
    positions = str(atoms.positions)

    return hash(symbols + positions)


def get_hash(string):
    return hashlib.md5(string.encode("utf-8")).hexdigest()


def get_center_of_mass(masses, positions):
    return masses @ positions / masses.sum()


def get_moments_of_inertia(positions, masses):
    """Get the moments of inertia along the principal axes.

    The three principal moments of inertia are computed from the
    eigenvalues of the symmetric inertial tensor. Periodic boundary
    conditions are ignored. Units of the moments of inertia are
    amu*angstrom**2.
    """
    com = get_center_of_mass(masses, positions)
    positions_ = positions - com  # translate center of mass to origin

    # Initialize elements of the inertial tensor
    I11 = I22 = I33 = I12 = I13 = I23 = 0.0
    for i in range(len(positions_)):
        x, y, z = positions_[i]
        m = masses[i]

        I11 += m * (y ** 2 + z ** 2)
        I22 += m * (x ** 2 + z ** 2)
        I33 += m * (x ** 2 + y ** 2)
        I12 += -m * x * y
        I13 += -m * x * z
        I23 += -m * y * z

    I = np.array([[I11, I12, I13], [I12, I22, I23], [I13, I23, I33]])

    evals, evecs = np.linalg.eigh(I)

    return evals
