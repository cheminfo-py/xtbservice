from rdkit import Chem
from rdkit.Chem import rdDistGeom
from ase import Atoms
from .cache import conformer_cache

def embed_conformer(mol, num_conformer: int = 10, prune_tresh: float = 0.1):
    """Use Riniker/Landrum conformer generator: https://pubs.acs.org/doi/10.1021/acs.jcim.5b00654"""
    param = rdDistGeom.ETKDGv2()
    param.pruneRmsThresh = prune_tresh
    cids = rdDistGeom.EmbedMultipleConfs(mol, num_conformer, param)
    # mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    # AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, mmffVariant='MMFF94s')


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
        mol = Chem.MolFromMolBlock(molfile, sanitize=False, removeHs=False)
        #Chem.calcImplicitValence(mol)
        mol.UpdatePropertyCache(strict=False)
        embed_conformer(mol)
        result =  rdkit2ase(mol)
        conformer_cache.set(molfile, result)
    return result

def smiles2ase(smiles: str) -> Atoms:
    try:
        result = conformer_cache.get(smiles)
    except KeyError:
        pass

    if result is None:
        mol = Chem.MolFromSmiles(smiles)
        refmol = Chem.AddHs(Chem.Mol(mol))
        embed_conformer(refmol)
        result =  rdkit2ase(refmol)
        conformer_cache.set(smiles, result)
    return result



def hash_atoms(atoms: Atoms) -> str: 
    symbols = str(atoms.symbols)
    positions = str(atoms.positions) 

    return hash(symbols + positions)