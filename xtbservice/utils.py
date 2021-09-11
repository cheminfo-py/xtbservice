from rdkit import Chem
from rdkit.Chem import rdDistGeom
from ase import Atoms

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


def smiles2ase(smiles: str) -> Atoms:
    mol = Chem.MolFromSmiles(smiles)
    refmol = Chem.AddHs(Chem.Mol(mol))
    embed_conformer(refmol)
    return rdkit2ase(refmol)
