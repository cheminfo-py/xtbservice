import hashlib

from ase import Atoms
from rdkit import Chem

from .cache import conformer_cache
from .conformers import embed_conformer


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
        mol = embed_conformer(mol)
        result =  rdkit2ase(mol), mol
        conformer_cache.set(molfile, result,  expire=None)
    return result

def smiles2ase(smiles: str) -> Atoms:
    try:
        result = conformer_cache.get(smiles)
    except KeyError:
        pass

    if result is None:
        mol = Chem.MolFromSmiles(smiles)
        refmol = Chem.AddHs(Chem.Mol(mol))
        refmol = embed_conformer(refmol)
        result =  rdkit2ase(refmol), refmol
        conformer_cache.set(smiles, result,  expire=None)
    return result



def hash_atoms(atoms: Atoms) -> str: 
    symbols = str(atoms.symbols)
    positions = str(atoms.positions) 

    return hash(symbols + positions)


def get_hash(string): 
    return hashlib.md5(string.encode("utf-8")).hexdigest()