from gemlog import *
import pytest
import gemlog
import gemlog.gemconvert
from gemlog.gem_network import summarize_gps_terminal
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

def test_make_inventory():
    gemlog.gemconvert.main(['-i', '../demo/demo/raw/'])
    gps_dir = 'gps'
    station_info_file = '../demo/demo/station_info.txt'
    summarize_gps_terminal([gps_dir, '-s', station_info_file])
