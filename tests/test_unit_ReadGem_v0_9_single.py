from gemlog import *
import pytest
## This function accounts for most of the runtime, so it's a candidate for speeding up.
## Any exception raised in this function should be interpreted as a bad data input raw file and
## handled immediately. Data files can be bad in unpredictable ways and can trigger exceptions
## in multiple places. It must never raise an exception when given a valid data file.

def test_ReadGem_v0_9_single():
    ## Check to make sure that the function does not crash when given good input, and that the
    ## output types and values are exactly correct. Flaws in the output (including hard-to-see
    ## flaws) cause hard-to-trace problems downstream.
    fn = 'data/FILE0001.077' 
    offset = 72000000.0 + 5263 # this is a reasonable offset to use for this file; can't use 0 because bugs might slip through
    eps = 1e-12
    L1 = ReadGem_v0_9_single(fn, offset)
    L0 = slow_ReadGem_v0_9_single(fn, offset)
    assert np.abs(L1['data'] - L0['data']).max() < eps
    for i in ['gps', 'metadata']:
        for x in L0[i].keys():
            assert np.abs(L0[i][x] - L1[i][x]).max() < eps
            assert L0[i][x].dtype == L1[i][x].dtype


def test_read_with_cython():
    ## Check to make sure that the function does not crash when given good input, and that the
    ## output types and values are exactly correct. Flaws in the output (including hard-to-see
    ## flaws) cause hard-to-trace problems downstream.
    fn = 'data/FILE0001.077' 
    offset = 72000000.0 + 5263 # this is a reasonable offset to use for this file; can't use 0 because bugs might slip through
    eps = 1e-5
    L1 = read_with_cython(fn, offset)
    L0 = slow_ReadGem_v0_9_single(fn, offset)
    assert np.abs(L1['data'] - L0['data']).max() < eps
    for i in ['gps', 'metadata']:
        for x in L0[i].keys():
            assert np.abs(L0[i][x] - L1[i][x]).max() < eps
            assert L0[i][x].dtype == L1[i][x].dtype


## These tests ensure that it raises an exception when reading bad raw
## files. No need to classify the type of exception; any exception
## should be interpreted as a bad file.
def test_ReadGem_v0_9_single_empty():
    with pytest.raises(EmptyRawFile):
        ReadGem_v0_9_single('data/FILE0000.000') # test empty file

def test_ReadGem_v0_9_single_corrupt():
    with pytest.raises(Exception):
        ReadGem_v0_9_single('data/FILE0023.096') # test a malformed file
