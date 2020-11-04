from gemlog.gemlog import _read_single_v0_8, EmptyRawFile, CorruptRawFileNoGPS, CorruptRawFile
from gemlog.gemlog import (
    _read_0_8_with_pandas, read_gem
)
import numpy as np
import pytest, shutil, os

import obspy

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

# use scope='session' to only evaluate this fixture once:
@pytest.fixture(scope='session')
def inputs():
    # this is a reasonable offset to use for this file;
    # can't use 0 because bugs might slip through
    offset = 72000000.0 + 5263
    return '../data/v0.8/raw/FILE0000.TXT', offset


@pytest.fixture(scope='session')
def test_read_single_v0_8(inputs):
    # serves as an implicit check that the reference reader doesn't error, but
    # would still be good to test its return values for correctness
    fn, offset = inputs
    return _read_single_v0_8(fn, offset)


def test_read_gem_integration():
    data_st = read_gem(path = '../data/v0.8/raw/', SN = '014')['data']
    data_st.merge()
    data = data_st[0]
    reference_st = obspy.read('../data/v0.8/converted_reference/*')
    reference_st.merge()
    reference = reference_st[0]

    ## make sure the waveforms span the same time window
    dt1 = data.stats.starttime
    rt1 = reference.stats.starttime
    dt2 = data.stats.endtime
    rt2 = reference.stats.endtime
    assert np.abs(dt1 - rt1) <= 0.01
    assert np.abs(dt2 - rt2) <= 0.01

    ## make sure the waveforms are similar. the methods are probably a bit different in the reference (from R.gemlog) so don't expect a perfect match
    t1 = np.max([dt1, rt1])
    t2 = np.min([dt2, rt2])
    eps = 1e-4
    data.trim(t1-eps, t2+eps)
    reference.trim(t1-eps, t2+eps)
    diff_ref = np.diff(reference.data)
    res_ref = np.min(np.abs(diff_ref[diff_ref != 0]))
    assert len(data.data) == len(reference.data)
    assert np.std(data.data - reference.data / res_ref) < 1
    
