import numpy as np
import pytest, shutil, os, obspy
import gemlog.gem_network

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

## example pytest decorators: see test_unit_read_functions for sample applications
# use scope='session' to only evaluate this fixture once (even if multiple tests call it)
#@pytest.fixture(scope='session')
# use parametrize to run the test with multiple inputs
#@pytest.mark.parametrize('reader_function', [_read_with_cython, _read_with_pandas])

def test_file_name_format():
    t1 = obspy.UTCDateTime('2020-01-03T05:02:08')
    outfile_pattern = '{year}-{mon}-{day}T{hour}_{min}_{sec}.{jd}.{net}.{sta}.{loc}.{chan}.{fmt}'
    trace_info_dict = {'year':t1.year, 'mon':t1.month, 'day':t1.day,
                       'jd':t1.julday, 'hour':t1.hour, 'min':t1.minute,
                       'sec':t1.second, 'net':'XP', 'sta':'PARK',
                       'loc':'00', 'chan':'HDF', 'fmt':'mseed'}
    outfile_pattern = gemlog.gem_network._fix_file_name_digits(outfile_pattern)
    output_file = outfile_pattern.format(**trace_info_dict)
    assert output_file == '2020-01-03T05_02_08.003.XP.PARK.00.HDF.mseed'
    return
