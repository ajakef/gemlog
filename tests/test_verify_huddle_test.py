import numpy as np
import pytest, shutil, os, obspy
import gemlog

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
    gemlog.huddle_test.verify_huddle_test('../data/test_data/huddle_test')
    # this test has disagreeing stop times because a battery died. add a check for that.
    # report is made in data folder. need to add an option to make the report in any folder.
    
