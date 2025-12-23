from gemlog.gps_debug import check_raw_file_for_gps_steps, find_gps_discontinuity
import numpy as np
import pytest

def test_check_raw_file_for_gps_steps():

    assert check_raw_file_for_gps_steps('/home/jake/Work/gemlog_python/data/test_data/gps_time_discontinuity/FILE0061.200') # has GPS step; should return True

    assert not check_raw_file_for_gps_steps('/home/jake/Work/gemlog_python/data/v1.10/FILE0002.232') # no GPS step; should return False

    #assert not check_raw_file_for_gps_steps('/home/jake/2024-04-26_RoofTestGPS/raw/FILE0004.356') # has a reversed GPS step; should return False


def test_find_gps_steps():
    assert np.round(find_gps_discontinuity('/home/jake/Work/gemlog_python/data/test_data/gps_time_discontinuity/FILE0061.200').reset_index().loc[0, 'msPPS'], 2) == 8867587.50 # previously used "index", but that's sensitive to how we filter out bad GPS points, which is not what we're testing here

    assert find_gps_discontinuity('/home/jake/Work/gemlog_python/data/v1.10/FILE0002.232').shape[0] == 0 # no steps
