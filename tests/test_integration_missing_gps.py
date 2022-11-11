import numpy as np
import pytest
import gemlog
import obspy
import sys
import os
import shutil
from gemlog import gem_cat, convert
from gemlog.core import *
from gemlog.core import _convert_one_file

def setup_module():
    try:
        shutil.rmtree('tmp') # might fail if the directory doesn't exist
    except:
        pass
    os.makedirs('tmp')
    os.makedirs('tmp/raw_merged')
    os.makedirs('tmp/converted')
    os.chdir('tmp')
    print(os.getcwd())
    
def teardown_module():
    os.chdir('..')
    shutil.rmtree('tmp')

def test_demo_gem_cat():
    ## following is drawn as directly as possible from demo/README.md
    gem_cat('../demo_missing_gps/raw_missing_gps/', './raw_merged', '077')
    convert(rawpath = './raw_merged', convertedpath = './converted', SN = '077')
    st = obspy.read('./converted/*')
    st_ref = obspy.read('../demo_missing_gps/converted_with_gps/*')
    st.merge()
    st_ref.merge()
    ## check the start and end times
    assert np.abs(st[0].stats.starttime - st_ref[0].stats.starttime) < 1
    assert np.abs(st[0].stats.endtime - st_ref[0].stats.endtime) < 1
    ## check the number of points
    assert np.abs(st[0].stats.npts - st_ref[0].stats.npts) < 10
    ## check the amplitudes
    assert np.abs(np.std(st[0].data)/np.std(st_ref[0].data) - 1) < 0.05
    
#test_demo_missing_gps()


## check that gemconvert_single handles potential file problems correctly
def test_demo_gemconvert_single():
    with pytest.raises(MissingRawFiles):
        _convert_one_file('../demo_missing_gps/raw_missing_gps/FILE9999.999') # does not exist

    with pytest.raises(CorruptRawFileNoGPS):
        _convert_one_file('../demo_missing_gps/raw_missing_gps/FILE0001.077', require_gps = True) # lacks GPS data

    _convert_one_file('../demo_missing_gps/raw_missing_gps/FILE0001.077', require_gps = False) # lacks GPS data
    _convert_one_file('../demo_missing_gps/raw_missing_gps/FILE0000.077', require_gps = True) # has GPS data


