from .conformer_generator import ConformerGenerator
from rdkit import Chem
from .models import Conformer, ConformerLibrary


def embed_conformer(
    mol, forcefield: str = "uff", num_conformer: int = 10, prune_tresh: float = 0.1
):
    """Use Riniker/Landrum conformer generator: https://pubs.acs.org/doi/10.1021/acs.jcim.5b00654"""
    conf_generator = ConformerGenerator()
    mol, _ = conf_generator.generate_conformers(mol)
    return mol


def conformers_from_smiles(smiles, forcefield, rmsd_threshold, max_conformers):
    mol = Chem.MolFromSmiles(smiles)
    return generate_conformers_from_mol(mol, forcefield, rmsd_threshold, max_conformers)


def conformers_from_molfile(molfile, forcefield, rmsd_threshold, max_conformers):
    mol = Chem.MolFromMolBlock(molfile, sanitize=True, removeHs=False)
    mol.UpdatePropertyCache(strict=False)

    return generate_conformers_from_mol(mol, forcefield, rmsd_threshold, max_conformers)


def generate_conformers_from_mol(mol, forcefield, rmsd_threshold, max_conformers):
    conf_generator = ConformerGenerator(
        max_conformers=max_conformers,
        force_field=forcefield,
        rmsd_threshold=rmsd_threshold,
    )
    mol, energies = conf_generator.generate_conformers(mol)
    conformers = []
    print(mol.GetNumConformers())
    for i in range(mol.GetNumConformers()):
        conf = mol.GetConformer(i)
        energy = energies[i]
        conformers.append(Conformer(molFile=Chem.MolToMolBlock(mol, confId=i), energy=energy))

    return ConformerLibrary(conformers=conformers)
