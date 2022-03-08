# -*- coding: utf-8 -*-
from cgitb import reset
import os

import numpy as np
from fastapi.testclient import TestClient

from xtbservice import __version__, app

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
client = TestClient(app)

from .help import clear_caches

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
    print(response.text)
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


def collect_node_types(d):
    counts = {"translation": 0, "rotation": 0, "vibration": 0}

    for n in range(len(d["modes"])):
        mode_type = d["modes"][n]["modeType"]

        if mode_type == "translation":
            counts["translation"] += 1
        elif mode_type == "rotation":
            counts["rotation"] += 1
        elif mode_type == "vibration":
            counts["vibration"] += 1

    return counts


def test_mode_type_assignent():
    # original bug report on dimethylacetamide
    clear_caches()
    response = client.get("/v1/ir?smiles=CC(N(C)C)=O")
    d = response.json()
    mode_types = collect_node_types(d)
    assert mode_types["translation"] == 3
    assert mode_types["rotation"] == 3

    # linear molecule
    clear_caches()
    response = client.get("/v1/ir?smiles=O=C=O")
    d = response.json()
    co2_mode_types = collect_node_types(d)
    assert co2_mode_types["translation"] == 3
    assert co2_mode_types["rotation"] == 2

    # water molecule
    clear_caches()
    response = client.get("/v1/ir?smiles=O")
    d = response.json()
    water_mode_types = collect_node_types(d)
    assert water_mode_types["translation"] == 3
    assert water_mode_types["rotation"] == 3


def test_raman():
    clear_caches()
    response = client.get("/v1/ir?smiles=O=C=O")
    d = response.json()
    assert np.abs(np.argmax(d["ramanIntensities"]) - 3224) < 5


def test_failing_raman():
    clear_caches()

    response = client.get('/v1/ir?smiles=Cl[Ti](Cl)(Cl)Cl')
    d = response.json()
    assert d['ramanIntensities'] == None