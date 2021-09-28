# -*- coding: utf-8 -*-
import os

from fastapi.testclient import TestClient

from xtbservice import __version__, app

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
client = TestClient(app)

from xtbservice.cache import (
    conformer_cache,
    ir_cache,
    ir_from_molfile_cache,
    ir_from_smiles_cache,
    opt_cache,
)


def clear_caches():
    for cache in [
        ir_from_smiles_cache,
        ir_from_molfile_cache,
        opt_cache,
        ir_cache,
        conformer_cache,
    ]:
        cache.clear()


def test_from_smiles():
    clear_caches()
    response = client.get("/v1/ir?smiles=CCCC")
    assert response.status_code == 200
    d = response.json()
    assert len(d["modes"]) == 3 * 16 - 6


def test_from_molfile():
    clear_caches()
    with open(os.path.join(THIS_DIR, "test_files", "molfile.mol"), "r") as handle:
        mol = handle.read()
    response = client.post("/v1/ir", json={"molFile": mol})
    assert response.status_code == 200
    d = response.json()
    assert 8 in d["modes"][35]["mostDisplacedAtoms"]
    assert 6 in d["modes"][35]["mostDisplacedAtoms"]
    assert 7 in d["modes"][35]["mostDisplacedAtoms"]


def test_most_contributing_atoms_bonds():
    clear_caches()
    with open(os.path.join(THIS_DIR, "test_files", "benzaldehyde.mol"), "r") as handle:
        mol = handle.read()
    response = client.post("/v1/ir", json={"molFile": mol})
    assert response.status_code == 200
    d = response.json()
    assert 1 in d["modes"][35]["mostContributingAtoms"]
    assert 2 in d["modes"][35]["mostContributingAtoms"]
    assert 8 in d["modes"][35]["mostContributingAtoms"]
    assert len(d["modes"][35]["mostContributingAtoms"]) == 3
