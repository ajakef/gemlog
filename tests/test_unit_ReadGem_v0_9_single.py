from gemlog import ReadGem_v0_9_single, EmptyRawFile
from gemlog.gemlog import (
    _read_with_cython, _read_with_pandas, _slow_ReadGem_v0_9_single
)
import numpy as np
import pytest
# Any exception raised in this function should be interpreted as a bad data
# input raw file and handled immediately. Data files can be bad in
# unpredictable ways and can trigger exceptions in multiple places. It must
# never raise an exception when given a valid data file.


@pytest.fixture(scope='session')
def inputs():
    # this is a reasonable offset to use for this file;
    # can't use 0 because bugs might slip through
    offset = 72000000.0 + 5263
    return 'data/FILE0001.077', offset


# use scope='session' to only evaluate this fixture once:
@pytest.fixture(scope='session')
def reference_output(inputs):
    # serves as an implicit check that the reference reader doesn't error, but
    # would still be good to test its return values for correctness
    fn, offset = inputs
    return _slow_ReadGem_v0_9_single(fn, offset)


def assert_gem_results_equal(L0, L1, eps=1e-12):
    assert np.abs(L1['data'] - L0['data']).max() < eps
    for i in ['gps', 'metadata']:
        for x in L0[i].keys():
            assert np.abs(L0[i][x] - L1[i][x]).max() < eps
            assert L0[i][x].dtype == L1[i][x].dtype


@pytest.mark.parametrize('reader_function', [_read_with_cython,
                                             _read_with_pandas,
                                             ReadGem_v0_9_single])
def test_file_readers(reader_function, inputs, reference_output):
    # test that the file readers work on a legitimate file,
    # and test that the output matches the reference output
    fn, offset = inputs
    actual_output = reader_function(fn, offset)
    assert_gem_results_equal(reference_output, actual_output)


# These tests ensure that it raises an exception when reading bad raw
# files. No need to classify the type of exception; any exception
# should be interpreted as a bad file.
@pytest.mark.parametrize('reader_function', [_read_with_cython,
                                             _read_with_pandas,
                                             ReadGem_v0_9_single])
def test_ReadGem_v0_9_single_empty(reader_function):
    with pytest.raises(EmptyRawFile):
        reader_function('data/FILE0000.000')  # test empty file


@pytest.mark.parametrize('reader_function', [_read_with_cython,
                                             _read_with_pandas,
                                             ReadGem_v0_9_single])
def test_ReadGem_v0_9_single_corrupt(reader_function):
    with pytest.raises(Exception):
        reader_function('data/FILE0023.096')  # test a malformed file
