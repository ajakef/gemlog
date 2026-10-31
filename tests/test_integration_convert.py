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
import numpy as np

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
    #shutil.rmtree('tmp')

## test a large block of files so that the loop in gemconvert is definitely covered by tests
def test_gemconvert_v110():
    gemlog.convert(rawpath='../data/v1.10/', convertedpath = 'mseed', SN= '232', blockdays = 0.5)
    st = obspy.read('mseed/*232..HDF.mseed')
    assert len(st) == 4
    shutil.rmtree('mseed/')

def test_gemconvert_AspenCSV001():
    gemlog.convert(rawpath='../data/AspenCSV/', convertedpath = 'mseed', SN= '977', blockdays = 0.5, nums = np.arange(20, 30)) # there's an unrelated .977 file; exclude it with nums
    ## converted data must keep read_gem's sample times, including after the data gap (16:13:36 to 16:14:05)
    st = obspy.read('mseed/*977..HDF.mseed')
    assert len(st) == 2
    raw = read_gem('../data/AspenCSV/', SN = '977', nums = np.arange(20, 30))['data'].select(channel = 'HDF')
    assert len(raw) == 2
    for tr, tr_raw in zip(sorted(st, key = lambda tr: tr.stats.starttime), raw):
        offset = (tr.stats.starttime - tr_raw.stats.starttime) / tr.stats.delta # samples
        assert abs(offset - np.round(offset)) < 0.001
    assert obspy.UTCDateTime('2025-10-20T16:14:05.063') in [tr.stats.starttime for tr in st] # was moved to .065
    shutil.rmtree('mseed/')
def test_gemconvert_AspenCSV01_24hr():
    gemlog.convert(rawpath='../data/test_data/24_hour_aspen_files/', convertedpath = 'mseed', SN= '10025') 
    shutil.rmtree('mseed/')

## test a block of files including a long data gap--e.g., old data was left on the disk and then new data recorded afterward
## if this works, it will only make mseed files for the 8 days when we actually have data
## if it fails, there will be many interpolated mseed files through the data gap(mid-april through mid-may)
def test_gemconvert_long_data_gaps():
    gemlog.convert(rawpath='../data/test_data/long_data_gaps/', convertedpath = 'mseed', SN= '128', blockdays = 0.5)
    st = obspy.read('mseed/*..HDF.mseed')
    assert len(st) == 8
    assert len(st.slice(obspy.UTCDateTime('2023-05-01'), obspy.UTCDateTime('2023-05-20'))) == 0
    shutil.rmtree('mseed/')
    
