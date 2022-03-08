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

def filter_modes(modes):
    return [mode for mode in modes if mode["modeType"] =='vibration']