# -*- coding: utf-8 -*-
import os

IMAGINARY_FREQ_THRESHOLD = os.getenv("IMAGINARY_FREQ_THRESHOLD", 10)
MAX_ATOMS = os.getenv("MAX_ATOMS", 80)
TIMEOUT = os.getenv("TIMEOUT", 500)
