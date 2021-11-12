# -*- coding: utf-8 -*-
import os

from fastapi.logger import logger

IMAGINARY_FREQ_THRESHOLD = int(os.getenv("IMAGINARY_FREQ_THRESHOLD", 10))
MAX_ATOMS = int(os.getenv("MAX_ATOMS", 50))
TIMEOUT = int(os.getenv("TIMEOUT", 100))

logger.info(
    f"Settings: IMAGINARY_FREQ_THRESHOLD: {IMAGINARY_FREQ_THRESHOLD}, MAX_ATOMS: {MAX_ATOMS}, TIMEOUT: {TIMEOUT}"
)
