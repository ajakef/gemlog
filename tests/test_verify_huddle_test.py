import numpy as np
import pytest, shutil, os, obspy
import gemlog
from gemlog.huddle_test import verify_huddle_test

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


def test_verify_huddle_test():
    verify_huddle_test('../data/test_data/huddle_test')
    # this test has disagreeing stop times because a battery died. add a check for that.
    # report is made in data folder. need to add an option to make the report in any folder.

def test_verify_huddle_test_aspen():
    verify_huddle_test('../data/test_data/huddle_test_aspen')
    # this test data uses 5-digit serial numbers and the Aspen metadata format (V/mA/RH
    # instead of A2/A3/FIFO/overruns/stack). report is made in data folder.

