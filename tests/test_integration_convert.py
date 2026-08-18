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
import pandas as pd

from gemlog.core import _convert_incremental

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
    shutil.rmtree('mseed/')


def test_incremental_conversion_assigns_midnight_boundary_once(monkeypatch):
    start = obspy.UTCDateTime('2024-01-01T23:59:59Z')
    traces = [obspy.Trace(
        data=np.ones(400, dtype=np.int32),
        header={'starttime': start, 'sampling_rate': 200,
                'station': 'TEST', 'network': '', 'location': '',
                'channel': 'HDF'})]

    def fake_read_gem(*args, **kwargs):
        number = kwargs['nums'][0]
        return {
            'data': obspy.Stream([traces[number - 1]]),
            'metadata': pd.DataFrame({'millis': [number]}),
            'gps': pd.DataFrame(),
            'header': pd.DataFrame({'bitweight_Pa': [1.0]}),
            'problems': [],
        }

    written = []
    calls = []

    def capture_window(stream, start, end, *args):
        calls.append((start, end, len(stream)))
        written.extend(stream.slice(start, end - 1e-6, nearest_sample=False).split())

    monkeypatch.setattr(gemlog.core, 'read_gem', fake_read_gem)
    monkeypatch.setattr(gemlog.core, '_write_utc_window', capture_window)
    _convert_incremental(
        rawpath='.', convertedpath='mseed', metadatapath='metadata',
        metadatafile='', gpspath='gps', gpsfile='', t1=-np.inf, t2=np.inf,
        nums=np.array([1]), SN='977', bitweight=np.nan, units='Pa',
        time_adjustment=0, blockdays=1, file_length_hour=1,
        station='TEST', network='', location='', output_format='MSEED')

    assert calls == [
        (obspy.UTCDateTime('2024-01-01T23:00:00Z'),
         obspy.UTCDateTime('2024-01-02T00:00:00Z'), 1),
        (obspy.UTCDateTime('2024-01-02T00:00:00Z'),
         obspy.UTCDateTime('2024-01-02T01:00:00Z'), 1),
    ]
    assert sum(len(trace) for trace in written) == 400
    assert sorted(trace.stats.starttime for trace in written) == [
        obspy.UTCDateTime('2024-01-01T23:59:59Z'),
        obspy.UTCDateTime('2024-01-02T00:00:00Z'),
    ]

## test a block of files including a long data gap--e.g., old data was left on the disk and then new data recorded afterward
## if this works, it will only make mseed files for the 8 days when we actually have data
## if it fails, there will be many interpolated mseed files through the data gap(mid-april through mid-may)
def test_gemconvert_long_data_gaps():
    gemlog.convert(rawpath='../data/test_data/long_data_gaps/', convertedpath = 'mseed', SN= '128', blockdays = 0.5)
    st = obspy.read('mseed/*..HDF.mseed')
    assert len(st) == 8
    assert len(st.slice(obspy.UTCDateTime('2023-05-01'), obspy.UTCDateTime('2023-05-20'))) == 0
    shutil.rmtree('mseed/')
    
