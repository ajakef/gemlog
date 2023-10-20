import gemlog
from gemlog import * # functions and classes not starting with _
import pytest
import shutil, os

def setup_module():
    try:
        shutil.rmtree('tmp') # exception if the directory doesn't exist
    except:
        pass
    os.makedirs('tmp')
    os.chdir('tmp')
    
def teardown_module():
    os.chdir('..')
    shutil.rmtree('tmp') 

## Rarely, the gps's time base can change for unknown reasons, resulting in an integer-second time
## offset. Make sure that these discontinuities are detected in the raw file where they occur
def test_step_detection():
    fn = '../data/gps_time_discontinuity/FILE0061.200'
    L = gemlog.core._read_single(fn, offset = 0)
    with pytest.raises(gemlog.core.CorruptRawFileDiscontinuousGPS):
        header_info = gemlog.core._calculate_drift(L, fn, require_gps = True)
