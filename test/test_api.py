# -*- coding: utf-8 -*-
from fastapi.testclient import TestClient

from xtb-service import __version__, app

client = TestClient(app)