from gemlog.gemlog import *
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


def test_read_gem_missing():
    with pytest.raises(MissingRawFiles):
        read_gem(np.arange(10, 20), '../data', SN = '000') # test missing files

def test_read_gem_missing():
    with pytest.raises(MissingRawFiles):
        read_gem(np.arange(5), '../data', SN = '000') # test purely empty files

def test_read_gem_edge_cases():
    print(os.getcwd())
    with pytest.raises(CorruptRawFile):
        read_gem(np.array([23]), '../data', SN = '096') # test a malformed file

def test_read_gem_good_data():
    ## read_gem always reads files in one block, so no sense in testing 25 files
    read_gem(np.arange(3), '../demo_missing_gps/raw_with_gps', SN = '077') # test good data


## Convert tests: ensure that it doesn't crash, and that the output mseed file is identical to a reference
def test_Convert_good_data():
    convert('../demo_missing_gps/raw_with_gps', SN = '077', convertedpath = 'test_output_mseed')
    output = obspy.read('test_output_mseed/2020-04-24T220000..077..HDF.mseed')
    reference = obspy.read('../demo_missing_gps/converted_with_gps/2020-04-24T220000..077..HDF.mseed')
    assert output.__eq__(reference)
