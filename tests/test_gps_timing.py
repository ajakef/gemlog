import gemlog
from gemlog import * # functions and classes not starting with _
from gemlog.gps_timing import _get_block_stats, get_GPS_spline, _check_step_within_block, _check_step_between_blocks
from gemlog.core import _calculate_drift
import pytest
import shutil, os
import numpy as np
default_deg1 = 0.001024

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


## Rarely, the gps's time base can change for unknown reasons, resulting in an integer-second time
## offset. Make sure that these discontinuities are detected in the raw file where they occur.
## Also, make sure that test_get_block_stats does its job otherwise.
def test_get_block_stats():
    slope_in = 0.001024 * (1-5e-5)
    x = np.arange(20)/slope_in + 1e6
    y = x * slope_in + 1e9
    #x = np.concatenate([np.arange(20), np.arange(100, 120)])

    # This dataset has no steps or spikes; should run without problems
    slope, x_mean, y_mean = _get_block_stats(x, y, default_deg1)
    assert approx_equal(slope, slope_in)
    assert approx_equal(x_mean, np.mean(x))
    assert approx_equal(y_mean, np.mean(y))
    assert not _check_step_within_block(x, y, default_deg1)[0]

    # This dataset has 1 positive & 1 negative spike; should run without problems and get same results
    true_y = x * slope_in + 1e9
    y = true_y + 10 * (x==3) - 15* (x==7)
    slope, x_mean, y_mean = _get_block_stats(x, y, default_deg1)
    assert approx_equal(slope_in, slope_in)
    assert approx_equal(x_mean, np.mean(x))
    assert approx_equal(y_mean, np.mean(true_y))
    assert not _check_step_within_block(x, y, default_deg1)[0]

    # This dataset has a step and should raise an exception
    y = (x > x[6]) + x * slope_in + 1e9
    assert _check_step_within_block(x, y, default_deg1)[0]
    with pytest.raises(gemlog.core.CorruptRawFileDiscontinuousGPS):
        slope, x_mean, y_mean = _get_block_stats(x, y, default_deg1)

def test_get_GPS_spline():
    # calculate a GPS spline and check that it handles outliers correctly
    fn = '../data/test_data/outlier_removal/FILE0011.202'
    L = gemlog.core._read_single(fn, offset = 0)
    G = L['gps']
    spline = _calculate_drift(L, fn, require_gps = True)['drift_spline']

    # check that it correctly extrapolates values outside the input range correctly (doesn't clamp them or have weird edge effects)
    ms = 0; assert np.abs(spline(ms) - 1693668082.114023) < 1e-6

    # check normal points between two blocks
    ms = 3e6; assert np.abs(spline(ms) - (1693668082.1136 + 1.02401829e-03*ms)) < 1e-4
    ms = 4e6; assert np.abs(spline(ms) - (1693668082.1136 + 1.02401829e-03*ms)) < 1e-4

    # check that it correctly ignores variation within a block and only connects block means
    assert np.std(spline(G['msPPS'][60:80]) - G['msPPS'][60:80] * 0.001024018) < (0.1 * np.std(G['t'][60:80] - G['msPPS'][60:80] * 0.001024018 - 1693668082.114023))
    
    
def test_integration_step_detection():
    ## this has a step between GPS blocks and should raise an exception
    fn = '../data/test_data/gps_time_discontinuity/FILE0061.200' # step between gps cycles
    L = gemlog.core._read_single(fn, offset = 0)
    with pytest.raises(gemlog.core.CorruptRawFileDiscontinuousGPS):
        header_info = gemlog.core._calculate_drift(L, fn, require_gps = True)

    fn = '../data/test_data/gps_time_discontinuity/FILE0004.356' # step within a gps cycle
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
#def test_robust_regress_recursion():
#    fn = '../data/test_data/outlier_removal/FILE0011.202'
#    L = gemlog.core._read_single(fn, offset = 0)
#    header_info = gemlog.core._calculate_drift(L, fn, require_gps = True)
#    assert approx_equal(header_info['drift_deg0'], 1.69366808e+09)
#    assert approx_equal(header_info['drift_deg1'], 1.02401829e-03)
#    assert approx_equal(header_info['drift_deg2'], -2.73652985e-17)
#    assert approx_equal(header_info['drift_deg3'], 3.17910536e-24)
