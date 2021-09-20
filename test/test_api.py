# -*- coding: utf-8 -*-
from fastapi.testclient import TestClient

from xtbservice import __version__, app

client = TestClient(app)

def test_from_smiles(): 
    response = client.get("/ir?smiles=CCCC")
    assert response.status_code == 200
    d =  response.json() 
    assert len(d['modes']) == 3 * 16  -6

# def test_from_molfile():
#     ...