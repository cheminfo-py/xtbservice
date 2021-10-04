# -*- coding: utf-8 -*-
from diskcache import Cache
import os

CACHE_DIR = os.getenv('CACHEDIR', 'ircache')

# 2 ** 30 = 1 GB 
opt_cache = Cache(CACHE_DIR, size_limit=2 ** 30, disk_min_file_size=0)
ir_cache = Cache(CACHE_DIR, size_limit=2 ** 30, disk_min_file_size=0)
conformer_cache = Cache(CACHE_DIR, size_limit=2 ** 30, disk_min_file_size=0)
ir_from_smiles_cache = Cache(CACHE_DIR, size_limit=2 ** 30, disk_min_file_size=0)
ir_from_molfile_cache = Cache(CACHE_DIR, size_limit=2 ** 30, disk_min_file_size=0)


if __name__ == "__main__":
    opt_cache.clear()
    ir_cache.clear()
    conformer_cache.clear()
    ir_from_smiles_cache.clear()
    ir_from_molfile_cache.clear()
