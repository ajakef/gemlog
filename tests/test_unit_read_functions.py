from gemlog.gemlog import EmptyRawFile, CorruptRawFileNoGPS, CorruptRawFile
from gemlog.gemlog import (
    _read_0_8_with_pandas, _read_with_pandas, _read_with_cython, read_gem, _read_single, _slow__read_single_v0_9, _process_gemlog_data
)
import numpy as np
import pytest, shutil, os, obspy


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
    return _read_single(fn, offset, version = '0.8')


def test_read_gem_integration():
    data_st = read_gem(path = '../data/v0.8/raw/', SN = '014')['data']
    data_st.merge()
    data = data_st[0]
    reference_st = obspy.read('../data/v0.8/converted_reference/*')
    reference_st.merge()
    reference = reference_st[0]
    #reference.stats.starttime += 0.01 # correction for mseed originally written with time offset

    ## make sure the waveforms span the same time window
    dt1 = data.stats.starttime
    rt1 = reference.stats.starttime
    dt2 = data.stats.endtime
    rt2 = reference.stats.endtime
    #assert np.abs(dt1 - rt1) <= 0.01
    #assert np.abs(dt2 - rt2) <= 0.01

    ## make sure the waveforms are similar. the methods are probably a bit different in the reference (from R.gemlog) so don't expect a perfect match
    t1 = np.max([dt1, rt1])
    t2 = np.min([dt2, rt2])
    eps = 1e-4
    data.trim(t1-eps, t2+eps)
    reference.trim(t1-eps, t2+eps)
    diff_ref = np.diff(reference.data)
    res_ref = np.min(np.abs(diff_ref[diff_ref != 0]))
    #assert len(data.data) == len(reference.data)
    #assert np.std(data.data - reference.data / res_ref) < 1

    ## timekeeping has changed enough since the R-gemlog days that these assertions are no longer useful, so commented out.
    
# Any exception raised in this function should be interpreted as a bad data
# input raw file and handled immediately. Data files can be bad in
# unpredictable ways and can trigger exceptions in multiple places. It must
# never raise an exception when given a valid data file.


@pytest.fixture(scope='session')
def inputs():
    # this is a reasonable offset to use for this file;
    # can't use 0 because bugs might slip through
    offset = 72000000.0 + 5263
    return '../demo_missing_gps/raw_with_gps/FILE0001.077', offset


# use scope='session' to only evaluate this fixture once:
@pytest.fixture(scope='session')
def reference_output(inputs):
    # serves as an implicit check that the reference reader doesn't error, but
    # would still be good to test its return values for correctness
    fn, offset = inputs
    return _slow__read_single_v0_9(fn, offset)

## test a good format 0.91 file to make sure it at least doesn't crash
@pytest.mark.parametrize('reader_function', [_read_with_cython, _read_with_pandas, _slow__read_single_v0_9])
def test_good_data_no_crash_0_91(reader_function):
    # serves as an implicit check that the reference reader doesn't error, but
    # would still be good to test its return values for correctness
    return reader_function('../data/v0.91/FILE0040.059', 5787)

def assert_gem_results_equal(L0, L1, eps=1e-12):
    assert np.abs(L1['data'] - L0['data']).max() < eps
    for i in ['gps', 'metadata']:
        for x in L0[i].keys():
            assert np.abs(L0[i][x] - L1[i][x]).max() < eps
            assert L0[i][x].dtype == L1[i][x].dtype


@pytest.mark.parametrize('reader_function', [_read_with_cython, _read_with_pandas])
def test_file_readers(reader_function, inputs, reference_output):
    # test that the file readers work on a legitimate file,
    # and test that the output matches the reference output
    fn, offset = inputs
    actual_output = reader_function(fn, offset)
    assert_gem_results_equal(reference_output, _process_gemlog_data(actual_output))


# These tests ensure that it raises an exception when reading bad raw
# files. No need to classify the type of exception; any exception
# should be interpreted as a bad file.
@pytest.mark.parametrize('reader_function', [_read_with_cython,
                                             _read_with_pandas,
                                             _read_single])
def test__read_single_v0_9_empty(reader_function):
    with pytest.raises(EmptyRawFile):
        reader_function('../data/FILE0000.000')  # test empty file


@pytest.mark.parametrize('reader_function', [_read_with_cython,
                                             _read_with_pandas,
                                             _slow__read_single_v0_9])
def test__read_single_v0_9_corrupt(reader_function):
    with pytest.raises(Exception): # can't use CorruptRawFile: it's up to _read_single to make that call.
        reader_function('../data/FILE0023.096')  # test a malformed file

@pytest.mark.parametrize('reader_function', [_read_with_cython,
                                             _read_with_pandas,
                                             _slow__read_single_v0_9])
def test__read_single_v0_9_no_gps(reader_function):
    with pytest.raises(CorruptRawFileNoGPS):
        reader_function('../demo_missing_gps/raw_missing_gps/FILE0001.077')  # test a missing-gps file
