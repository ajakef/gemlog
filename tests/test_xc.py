from gemlog import *
import numpy as np
import pandas as pd
import pytest
import gemlog
import gemlog.xcorr
import obspy
from obspy.clients.nrl import NRL
import sys
import os
import shutil
from obspy.core.inventory import Network, Station, Channel

def setup_module():
    try:
        shutil.rmtree('tmp') # might fail if the directory doesn't exist
    except:
        pass
    os.makedirs('tmp')
    os.chdir('tmp')
    print(os.getcwd())
    
    # make_data for test:
    km2deg = 9e-3

    ## scenario with wave propagating straight NE, so it arrives at station 000 first, station 001 0.06 s later, and station 002 just barely after that (0.02 s)
    x = np.array([0, 0, 0.040])
    y = np.array([0, 0.030, 0])

    G = np.array([[x[0] - x[1], y[0] - y[1]], 
                  [x[1] - x[2], y[1] - y[2]]])

    s = np.array([2, 2]) # s/km: 

    lags = G @ s
    stations = ['000', '001', '002']
    t = np.arange(1100) * 0.01 - 5
    def f(t):
        return np.sin(2*np.pi*7.1*t) * np.sin(2*np.pi*2*t) + np.random.normal(0, 0.001, 1100)
    
    st = obspy.Stream([
        obspy.Trace(data = f(t)),
        obspy.Trace(data = f(t+lags[0])),
        obspy.Trace(data = f(t+lags[1]+lags[0]))])
    for i, tr in enumerate(st):
        tr.stats.delta = 0.01
        tr.stats.station = stations[i]
        tr.write(tr.stats.starttime.isoformat().replace(':', '_') + '.' + tr.id + '.mseed')
    st.write('tmp.mseed')
    
    inv = obspy.Inventory()
    inv.networks.append(Network('', stations = [
        Station(stations[0], latitude = 0,longitude = 0, elevation = 0, channels = [Channel('HDF', '', latitude = 0,longitude = 0, elevation = 0, depth = 0)]),
        Station(stations[1], latitude = y[1]*km2deg,longitude = x[1]*km2deg, elevation = 0, channels = [Channel('HDF', '', latitude = y[1]*km2deg,longitude = x[1]*km2deg, elevation = 0, depth = 0)]),
        Station(stations[2], latitude = y[2]*km2deg, longitude = x[2] * km2deg, elevation = 0, channels = [Channel('HDF', '', latitude = y[2]*km2deg,longitude = x[2]*km2deg, elevation = 0, depth = 0)])
        ]))

    inv.write('tmp_inv.xml', format='stationxml')

def teardown_module():
    os.chdir('..')
    shutil.rmtree('tmp')

def test_xcorr():
    mseed_path_pattern = '1970*'
    t1 = '19700101_000000'
    t2 = '19700101_000010'
    gemlog.xcorr.xcorr_all_terminal([mseed_path_pattern, '-1', t1, '-2', t2, '-o', 'xc_results.csv', '-u', '4'])
    results = pd.read_csv('xc_results.csv')
    assert all(np.abs(results['lag_.000._.001.'] + 0.06) < 0.01) 
    assert all(np.abs(results['lag_.001._.002.'] + 0.02) < 0.01) 
    assert all(np.abs(results['lag_.002._.000.'] - 0.08) < 0.01) 
    assert all(np.abs(results['lag_.000._.001.'] + results['lag_.001._.002.'] + results['lag_.002._.000.']) < 0.005) 

def test_calculate_direction():
    xcorr_df = pd.read_csv('xc_results.csv')

    gemlog.xcorr.calculate_direction_terminal([
        '-i', 'xc_results.csv', '-l', 'tmp_inv.xml', '-o', 'directions.csv'
    ])
    directions = pd.read_csv('directions.csv')

    assert all(np.abs(directions.backazimuth + 135) < 3) # back-azimuth should be within 3 degrees of correct -135
    assert all(np.abs(directions.slowness/np.sqrt(8) - 1) < 0.05) # slowness should be correct to within 5%
    
