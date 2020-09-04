from gemlog import *
import pytest

def test_ReadGem_missing():
    with pytest.raises(MissingRawFiles):
        ReadGem(np.arange(10, 20), 'data', SN = '000') # test missing files

def test_ReadGem_missing():
    with pytest.raises(MissingRawFiles):
        ReadGem(np.arange(5), 'data', SN = '000') # test purely empty files

def test_ReadGem_edge_cases():
    with pytest.raises(CorruptRawFile):
        ReadGem(np.array([23]), 'data', SN = '096') # test a malformed file

def test_ReadGem_good_data():
    ## ReadGem always reads files in one block, so no sense in testing 25 files
    ReadGem(np.arange(3), 'data', SN = '077') # test good data


## Convert tests: ensure that it doesn't crash, and that the output mseed file is identical to a reference
def test_Convert_good_data():
    Convert('data', SN = '077', convertedpath = 'test_output_mseed')
    output = obspy.read('test_output_mseed/2020-04-24T22:00:00..077..HDF.mseed')
    reference = obspy.read('data/2020-04-24T22:00:00..077..HDF.mseed')
    assert output.__eq__(reference)
