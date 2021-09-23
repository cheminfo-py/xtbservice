from xtbservice.conformers import conformers_from_smiles

def test_conformer_generation(): 
    conf_library = conformers_from_smiles('Fc1ccc(cc1)[C@@]3(OCc2cc(C#N)ccc23)CCCN(C)C', 'uff', 0.5, 10)
    assert len(conf_library.conformers) > 2 
    assert conf_library.conformers[0].energy < conf_library.conformers[1].energy