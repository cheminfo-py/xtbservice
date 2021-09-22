# -*- coding: utf-8 -*-
"""Test if the API responds in the Docker image"""
import sys

import requests


r = requests.get("http://localhost:8091/v1/app_version/")
keys = r.json().keys()
if "version" in keys:
    print("OK")
    sys.exit(0)
else:
    print("error")
    sys.exit(1)
