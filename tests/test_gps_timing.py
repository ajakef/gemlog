import gemlog
from gemlog import * # functions and classes not starting with _
import pytest
import shutil, os
import numpy as np

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

def approx_equal(x, y, d = 0.001):
    return np.abs((x - y)/x) < d


def test_get_block_stats():
    x = np.arange(10)
    # This dataset has no steps or spikes; should run without problems
    slope, x_mean, y_mean = _get_block_stats(x, y = x * 0.0001)
    assert approx_equal(slope, 0.0001)
    assert approx_equal(x_mean, 4.5)
    assert approx_equal(y_mean, 0.00045)

    # This dataset has 2 spikes; should run without problems and get same results
    slope, x_mean, y_mean = _get_block_stats(x, y = x * 0.0001 + 10 * ((x==3) | (x==7)))
    assert approx_equal(slope, 0.0001)
    assert approx_equal(x_mean, 4.5)
    assert approx_equal(y_mean, 0.00045)

    # This dataset has a step and should raise an exception
    with pytest.raises(gemlog.core.CorruptRawFileDiscontinuousGPS):
        slope, x_mean, y_mean = _get_block_stats(x, y = (x > 6) + x * 0.0001)


## Rarely, the gps's time base can change for unknown reasons, resulting in an integer-second time
## offset. Make sure that these discontinuities are detected in the raw file where they occur
def test_unit_step_detection():
    ## this has a step and should raise an exception
    fn = '../data/test_data/gps_time_discontinuity/FILE0061.200'

    ## test a file with a spike

    
def test_integration_step_detection():
    ## this has a step and should raise an exception
    fn = '../data/test_data/gps_time_discontinuity/FILE0061.200'
    L = gemlog.core._read_single(fn, offset = 0)
    with pytest.raises(gemlog.core.CorruptRawFileDiscontinuousGPS):
        header_info = gemlog.core._calculate_drift(L, fn, require_gps = True)

    ## this does not have a step and should run without error
    fn = '../data/v1.10/FILE0002.232'
    L = gemlog.core._read_single(fn, offset = 0)
    header_info = gemlog.core._calculate_drift(L, fn, require_gps = True)

    ## test a file with a spike that should be ignored
    
    ## this has a reversed GPS step that should be ignored and run without error
    #fn = '/home/jake/2024-04-26_RoofTestGPS/raw/FILE0004.356' ## update this
    #L = gemlog.core._read_single(fn, offset = 0)
    #header_info = gemlog.core._calculate_drift(L, fn, require_gps = True)    
    
## _robust_regress works by recursively removing outliers from a cubic fit until there are none with z>4 or max abs dev > 0.01. Here's an example file with outliers where recursion is required; make sure it gives the right answer and doesn't break.
def test_robust_regress_recursion():
    fn = '../data/test_data/outlier_removal/FILE0011.202'
    L = gemlog.core._read_single(fn, offset = 0)
    header_info = gemlog.core._calculate_drift(L, fn, require_gps = True)
    assert approx_equal(header_info['drift_deg0'], 1.69366808e+09)
    assert approx_equal(header_info['drift_deg1'], 1.02401829e-03)
    assert approx_equal(header_info['drift_deg2'], -2.73652985e-17)
    assert approx_equal(header_info['drift_deg3'], 3.17910536e-24)
