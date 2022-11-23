## Right now, this test just makes sure that the demo runs without errors. In the future, it could be expanded to check the values of the test outputs--a much stronger test.

from gemlog import *
import pytest
import gemlog
import gemlog.gemconvert
import obspy
from obspy.clients.nrl import NRL
import sys
import os
import shutil

def setup_module():
    try:
        shutil.rmtree('tmp') # might fail if the directory doesn't exist
    except:
        pass
    os.makedirs('tmp')
    os.chdir('tmp')
    print(os.getcwd())
    
def teardown_module():
    os.chdir('..')
    shutil.rmtree('tmp')

## test a large block of files so that the loop in gemconvert is definitely covered by tests
def test_gemconvert_v110():
    gemlog.convert(rawpath='../data/v1.10/', convertedpath = 'mseed', SN= '232')
    st = obspy.read('mseed/*232..HDF.mseed')
    assert len(st) == 4
