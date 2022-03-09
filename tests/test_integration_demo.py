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

@pytest.mark.filterwarnings('ignore')
def test_demo():
    ## following is drawn as directly as possible from demo/README.md
    gemlog.gemconvert.main(['-i', '../demo/demo/raw/'])
    ##########################
    coords = gemlog.summarize_gps('gps', output_file = 'project_coords.csv', station_info = '../demo/demo/station_info.txt')
    gemlog.rename_files('mseed/*', station_info = '../demo/demo/station_info.txt', output_dir = 'renamed_mseed')
    #nrl = NRL()
    #response = nrl.get_response(sensor_keys = ['Gem', 'Gem Infrasound Sensor v1.0'],
    #datalogger_keys = ['Gem', 'Gem Infrasound Logger v1.0',
    #'0 - 128000 counts/V']) # may cause warning--ok to ignore
    response = gemlog.get_gem_response()
    ## manually add elevation to 'coords'. raise an issue on github if you know an
    ## easy-to-install, cross-platform way to automate this!
    coords['elevation'] = [1983, 1983, 1988, 1983, 1986, 1987] # from google earth
    
    ## create an inventory of all sensors used in this project--may cause warnings
    inv = gemlog.make_gem_inventory('../demo/demo/station_info.txt', coords, response)
    inv.write('NM_inventory.xml', format='STATIONXML')
    
    ## read the data
    data = obspy.read('renamed_mseed/*')
    print(data)

    ## combine traces so that each station has one trace
    data.merge()
    print(data)

    ## deconvolve the instrument responses using the inventory already created
    inv = obspy.read_inventory('NM_inventory.xml')
    data.remove_response(inv) # may cause warnings--ok to ignore
    
    ## filter data above 1 Hz (lower frequencies are often wind noise)
    data.filter("highpass", freq=1.0)
    
    ## trim the data around a known event
    t1 = obspy.UTCDateTime('2020-05-10T12:14:00')
    t2 = obspy.UTCDateTime('2020-05-10T12:15:00')
    data.trim(t1, t2)
    
    ## plot the results
    #data.plot() # suppress plotting

