from xtbservice.ir import ir_from_smiles
from .help import clear_caches, filter_modes
import pytest 

def get_ir_modes(modes, threshold=0.01):
    return [mode for mode in modes if mode['intensity'] > threshold]

def get_raman_modes(modes, threshold=0.01):
    return [mode for mode in modes if mode['ramanIntensity'] > threshold]

@pytest.mark.parametrize("technique", ["GFN2xTB", "GFNFF"])
def test_co2(technique):
    """3N-5 with center of inversion"""
    clear_caches()
    ir = ir_from_smiles('C(=O)=O', technique)

    modes = filter_modes(ir.modes)

    raman_intensities = get_raman_modes(modes)
    assert len(raman_intensities) == 1

    ir_intensities = get_ir_modes(modes)
    assert len(ir_intensities) == 3



@pytest.mark.parametrize("technique", ["GFN2xTB", "GFNFF"])
def test_cos(technique):
    """3N-5 no center of symmetry"""
    print(technique)
    clear_caches()
    ir = ir_from_smiles('C(=O)=S', technique)

    modes = filter_modes(ir.modes)

    raman_intensities = get_raman_modes(modes)
    assert len(raman_intensities) == 4

    ir_intensities =get_ir_modes(modes)
    assert len(ir_intensities) == 4


@pytest.mark.parametrize("technique", ["GFN2xTB", "GFNFF"])
def test_h2o(technique):
    """3N-6 DOF"""
    clear_caches()
    ir = ir_from_smiles('O', technique)

    modes = filter_modes(ir.modes)

    raman_intensities = get_raman_modes(modes)
    assert len(raman_intensities) == 3

    ir_intensities =get_ir_modes(modes)
    assert len(ir_intensities) == 3