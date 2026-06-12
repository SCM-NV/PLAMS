import pytest
import numpy as np
import shutil
from pathlib import Path 

from scm.plams.tools.converters import  vasp_output_to_ams
from scm.plams.tools.kftools import KFReader#, KFFile, KFHistory, KFTypedReadError

@pytest.fixture
def single_point_vasp(vasp_folder):
    return vasp_folder / "SPC"

@pytest.fixture
def npt_vasp(vasp_folder):
    return vasp_folder / "NpT"

@pytest.fixture
def relaxation_vasp(vasp_folder):
    return vasp_folder / "relaxation"

class TestVaspConverter:
    "Test Vasp Converter"  

    def test_singlepoint_conversion(self, single_point_vasp):
        rkf_dir = Path(vasp_output_to_ams(str(single_point_vasp)))
        reader = KFReader(rkf_dir/"ams.rkf")

        assert reader.read("General", "task") == "singlepoint" 

        shutil.rmtree(rkf_dir)
    
    def test_NpT_conversion(self, npt_vasp):
        rkf_dir = Path(vasp_output_to_ams(str(npt_vasp)))
        reader = KFReader(rkf_dir/"ams.rkf")

        assert reader.read("General", "task") == "moleculardynamics" 
        assert reader.read("MDHistory", "Time(1)") == list(np.arange(0, 10, dtype=float))
        
        shutil.rmtree(rkf_dir)

    @pytest.mark.skip(reason="ase can not recognize task; later version might")
    def test_relaxation_conversion(self, relaxation_vasp):
        rkf_dir = Path(vasp_output_to_ams(str(relaxation_vasp)))
        reader = KFReader(rkf_dir/"ams.rkf")

        assert reader.read("General", "task") == "geometryoptimization" 

        shutil.rmtree(rkf_dir)
    


    

