# -*- coding: utf-8 -*-
import os

from fastapi.logger import logger

IMAGINARY_FREQ_THRESHOLD = int(os.getenv("IMAGINARY_FREQ_THRESHOLD", 10))
MAX_ATOMS_XTB = int(os.getenv("MAX_ATOMS_XTB", 60))
MAX_ATOMS_FF = int(os.getenv("MAX_ATOMS_FF", 100))
TIMEOUT = int(os.getenv("TIMEOUT", 100))

logger.info(
    f"Settings: IMAGINARY_FREQ_THRESHOLD: {IMAGINARY_FREQ_THRESHOLD}, MAX_ATOMS_XTB: {MAX_ATOMS_XTB}, MAX_ATOMS_FF: {MAX_ATOMS_FF}, TIMEOUT: {TIMEOUT}"
)
