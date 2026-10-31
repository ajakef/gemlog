from gemlog.core import _apply_segments, _interp_time, _read_several, _merge_gaps
from gemlog.core import * # * doesn't load functions that start with _
import pytest
import shutil, os

import platform
import sys
import numpy as np
import scipy

print(sys.version)
print(platform.platform())
print(np.__version__)
print(scipy.__version__)

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

## unit test for key merging function in gemlog.core
def test_unit_merge_gaps():
    tr_1 = obspy.read()[0]
    tr_2 = obspy.read()[0]
    tr_2.stats.starttime = tr_1.stats.endtime+0.01
    result = _merge_gaps(obspy.Stream([tr_1, tr_2]))
    assert len(result) == 1
    tr_2.stats.starttime = tr_1.stats.endtime+0.03
    result = _merge_gaps(obspy.Stream([tr_1, tr_2]))
    assert len(result) == 1
    # this test is excessive for now; allowing wider gaps is probably fine
    tr_2.stats.starttime = tr_1.stats.endtime+0.04
    result = _merge_gaps(obspy.Stream([tr_1, tr_2]))
    assert len(result) == 2

## data after a gap must keep its own times, not be moved onto the sample grid of the data before the gap
def test_unit_merge_gaps_keeps_times():
    tr_1 = obspy.read()[0] # 100 sps
    for gap in [0.0632, 0.0681, 0.5047, 86400*30 + 0.0047]: # off tr_1's sample grid; the last is a month-long gap
        tr_2 = obspy.read()[0]
        tr_2.stats.starttime = tr_1.stats.endtime + gap
        result = _merge_gaps(obspy.Stream([tr_1.copy(), tr_2.copy()]))
        assert len(result) == 2
        assert result[0].stats.starttime == tr_1.stats.starttime
        assert result[1].stats.starttime == tr_2.stats.starttime
        assert result[0].stats.npts == tr_1.stats.npts
        assert result[1].stats.npts == tr_2.stats.npts

    ## small gaps are still interpolated, and the segment after a big gap keeps its times
    tr_2 = obspy.read()[0]
    tr_2.stats.starttime = tr_1.stats.endtime + 0.02
    tr_3 = obspy.read()[0]
    tr_3.stats.starttime = tr_2.stats.endtime + 0.03
    tr_4 = obspy.read()[0]
    tr_4.stats.starttime = tr_3.stats.endtime + 1.0632
    result = _merge_gaps(obspy.Stream([tr_1, tr_2, tr_3, tr_4]))
    assert len(result) == 2
    assert result[0].stats.starttime == tr_1.stats.starttime
    assert result[0].stats.npts == 3 * tr_1.stats.npts + 1 + 2 # 1 and 2 interpolated samples
    assert result[1].stats.starttime == tr_4.stats.starttime


## This test makes sure that spurious time gaps do not show up in converted data due to small
## timing errors in the conversion. The issue may have to do with max_step/min_step in
## gemlog.core._interp_time (which can make the function trigger-happy in finding data gaps) or
## in gemlog.core._robust_regress (which calculates a cubic regression between GPS-millis times).
## Spurious data gaps appear between files due to discrepancies between their drift corrections.

## All pairs of raw files added to data/spurious_data_gaps should have resulted in a spurious data
## gap before. Files should always be in pairs.

## if a file has weird timing, try code like this (run from /home/jake/Work/gemlog_python/data/test_data/spurious_data_gaps/2022-01-19_169_4-5)
# gemlog.core._plot_drift(sorted(os.listdir('.')))

def test_all_spurious_data_gaps():
    print(os.listdir('../data/test_data'))
    data_root = '../data/test_data/spurious_data_gaps/'
    dirs = sorted(os.listdir(data_root))
    for testdir in [os.path.join(data_root, d) for d in dirs]:
        file_list = sorted(os.listdir(testdir))
        L = _read_several([os.path.join(testdir, f) for f in file_list])
        D = L['data']
        D = np.hstack((D, _apply_segments(D[:,0], L['header']).reshape([D.shape[0],1])))
        st = _interp_time(D)
        if len(st) > 1:
            raise Exception('spurious data gap in ' + testdir)
        elif len(st) == 0:
            raise Exception('failed to convert data in ' + testdir)
        else:
            print(testdir + ' ok')
    
