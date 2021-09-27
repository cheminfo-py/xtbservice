# -*- coding: utf-8 -*-
"""
xtb-service
webservice providing xtb calculations
"""
# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

# Add imports here
from .xtbservice import *
