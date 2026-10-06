from gemlog.core import EmptyRawFile, CorruptRawFileNoGPS, CorruptRawFile
from gemlog.parsers import parse_gemfile
from gemlog.core import (
    _read_0_8_pd, _read_with_pandas, _read_with_cython, read_gem, _read_single,
    _read_several, _slow__read_single_v0_9, _process_gemlog_data, _read_SN,
    _read_format_version, _read_config_gem, _read_config_aspen
)
import numpy as np
import pytest, shutil, os, obspy


def setup_module():
    try:
        shutil.rmtree('tmp') # exception if the directory doesn't exist
    except:
        pass
    os.makedirs('tmp')
    os.chdir('tmp')
    
def teardown_module():
    os.chdir('..')
    shutil.rmtree('tmp', ignore_errors = True) 

# use scope='session' to only evaluate this fixture once:
@pytest.fixture(scope='session')
def inputs():
    # this is a reasonable offset to use for this file;
    # can't use 0 because bugs might slip through
    offset = 72000000.0 + 5263
    return '../data/v0.8/raw/FILE0000.TXT', offset

def test_parser():
    x = parse_gemfile(b'../data/v1.10/FILE0001.210')
    assert x[0][0,0] == 635
    assert x[2][0] == 7174
    x = parse_gemfile(b'../data/AspenCSV/FILE0009.004', n_channels = 4) # very basic aspen file, no compression, bad millis rollover. Don't try to convert it.
    assert x[0][0,0] == 49242
    assert x[2][0] == 6690
    assert np.sum(x[1] == b'D') == 54400
    x = parse_gemfile(b'../data/AspenCSV/FILE0001.017', n_channels = 4) # 2026-07-22 4-channel aspen file with compression. 

    # aspen file of format 0.1 with normal compression, with multiple compression formats possible
    x, linetype, millis = parse_gemfile(b'../data/AspenCSV/FILE0077.977', n_channels = 1, dt_ms = 5)
    w = np.where(linetype == b'D')[0]
    data = x[w,0]
    # lines being parsed shown in comments
    assert millis[w[6670]] == 2393 # 13
    assert data[6670] == 13
    assert millis[w[6671]] == 2399 # zk
    assert data[6671] == -2
    assert millis[w[6672]] == 2409 # D2409,q
    assert data[6672] == 4
    assert millis[w[6673]] == 2414 # D2409,q
    assert data[6673] == 4

    # 24-hour aspen file of format 0.1
    x, linetype, millis = parse_gemfile(b'../data/test_data/24_hour_aspen_files/FILE0001.025', n_row = int(200*3600*24*1.1), n_channels = 1, dt_ms = 5)
    
    
    

    #x = parse_gemfile(b'../data/AspenCSV/FILE0108.006') # more recent 4-channel test file with compression (ERB backyard)
    #assert x[0][0,0] == -5868
    #assert x[2][0] == 1015
    #assert np.sum(x[1] == b'D') == 136564
    
def test_read_SN():
    _read_SN('../data/v0.91/FILE0040.059')
    _read_SN('../data/v1.10/FILE0001.210')
    _read_SN('../data/AspenCSV/FILE0020.977') # very basic aspen file

def test_read_format():
    _read_format_version('../data/v0.91/FILE0040.059')
    _read_format_version('../data/v1.10/FILE0001.210')
    _read_format_version('../data/AspenCSV/FILE0020.977') # very basic aspen file

def test_read_config():
    _read_config_gem('../data/v0.91/FILE0040.059')
    _read_config_gem('../data/v1.10/FILE0001.210')
    _read_config_aspen('../data/AspenCSV/FILE0020.977') # very basic aspen file--replace with one that has a proper config line
    x = _read_config_aspen('../data/aspen_config/FILE0078.977')
    assert x['n_channels'] == 1
    assert x['sample_int_ms'] == 5
    x = _read_config_aspen('../data/aspen_config/FILE0079.977')
    assert x['n_channels'] == 1
    assert x['sample_int_ms'] == 5
    x = _read_config_aspen('../data/aspen_config/FILE0080.977')
    assert x['n_channels'] == 1
    assert x['sample_int_ms'] == 5
@pytest.fixture(scope='session')

def test_read_single_v0_8(inputs):
    # serves as an implicit check that the reference reader doesn't error, but
    # would still be good to test its return values for correctness
    fn, offset = inputs
    _read_single(fn, offset, version = '0.8')

def test_read_single_v0_8_NMEA_single_character_lines():
    _read_0_8_pd('../data/v0.8_NMEA_single_character_lines/raw/FILE0014.TXT')

def test_read_gem_integration():
    data_st = read_gem(path = '../data/v0.8/raw/', SN = '014')['data']
    data_st.merge()
    data = data_st[0]
    reference_st = obspy.read('../data/v0.8/converted_reference/*')
    reference_st.merge()
    reference = reference_st[0]
    #reference.stats.starttime += 0.01 # correction for mseed originally written with time offset

    ## make sure the waveforms span the same time window
    dt1 = data.stats.starttime
    rt1 = reference.stats.starttime
    dt2 = data.stats.endtime
    rt2 = reference.stats.endtime
    #assert np.abs(dt1 - rt1) <= 0.01
    #assert np.abs(dt2 - rt2) <= 0.01

    ## make sure the waveforms are similar. the methods are probably a bit different in the reference (from R.gemlog) so don't expect a perfect match
    t1 = np.max([dt1, rt1])
    t2 = np.min([dt2, rt2])
    eps = 1e-4
    data.trim(t1-eps, t2+eps)
    reference.trim(t1-eps, t2+eps)
    diff_ref = np.diff(reference.data)
    res_ref = np.min(np.abs(diff_ref[diff_ref != 0]))
    #assert len(data.data) == len(reference.data)
    #assert np.std(data.data - reference.data / res_ref) < 1

    ## timekeeping has changed enough since the R-gemlog days that these assertions are no longer useful, so commented out.
    
# Any exception raised in this function should be interpreted as a bad data
# input raw file and handled immediately. Data files can be bad in
# unpredictable ways and can trigger exceptions in multiple places. It must
# never raise an exception when given a valid data file.

def test_read_single_aspen():
    L = _read_single('../data/AspenCSV/FILE0020.977', version = 'AspenCSV0.01')
    #L = _read_single('../data/AspenCSV/FILE0108.006', version = 'AspenCSV0.01')
    #assert np.shape(L['data'])[0] == 136564
    x = _read_single('../data/AspenCSV/FILE0001.017', version = 'AspenCSV0.1') # 4-channel, recorded in lab in 2026-07
    x = _read_single('../data/test_data/24_hour_aspen_files/FILE0001.025', version = 'AspenCSV0.1')


def test_unit_read_several_aspen():
    L = _read_several(['../data/AspenCSV/FILE0020.977'], version = 'AspenCSV0.01')
    #L = _read_several(['../data/AspenCSV/FILE0108.006'], version = 'AspenCSV0.01')
    #assert np.shape(L['data'])[0] == 136564
    x = _read_several(['../data/AspenCSV/FILE0001.017'], version = 'AspenCSV0.1') # 4-channel, recorded in lab in 2026-07

    x = _read_several(['../data/test_data/24_hour_aspen_files/FILE0001.025'], version = 'AspenCSV0.1')
    x = _read_several(['../data/test_data/24_hour_aspen_files/FILE0001.025', '../data/test_data/24_hour_aspen_files/FILE0002.025'], version = 'AspenCSV0.1')
    

def test_read_gem_aspen():
    data_st = read_gem(path = '../data/AspenCSV/', nums = [1], SN = '10017')['data']
    data_st = read_gem(path = '../data/test_data/24_hour_aspen_files/', SN = '10025')['data']

@pytest.fixture(scope='session')
def inputs():
    # this is a reasonable offset to use for this file;
    # can't use 0 because bugs might slip through
    offset = 72000000.0 + 5263
    return '../demo_missing_gps/raw_with_gps/FILE0001.077', offset


# use scope='session' to only evaluate this fixture once:
@pytest.fixture(scope='session')
def reference_output(inputs):
    # serves as an implicit check that the reference reader doesn't error, but
    # would still be good to test its return values for correctness
    fn, offset = inputs
    return _slow__read_single_v0_9(fn, offset)

## test a good format 0.91 file to make sure it at least doesn't crash
@pytest.mark.parametrize('reader_function', [_read_with_cython, _read_with_pandas, _slow__read_single_v0_9])

def test_good_data_no_crash_0_91(reader_function):
    # serves as an implicit check that the reference reader doesn't error, but
    # would still be good to test its return values for correctness
    reader_function('../data/v0.91/FILE0040.059', 5787)

def assert_gem_results_equal(L0, L1, eps=1e-12):
    assert np.abs(L1['data'] - L0['data']).max() < eps
    for i in ['gps', 'metadata']:
        for x in L0[i].keys():
            assert np.abs(L0[i][x] - L1[i][x]).max() < eps
            assert L0[i][x].dtype == L1[i][x].dtype


@pytest.mark.parametrize('reader_function', [_read_with_cython, _read_with_pandas])
def test_file_readers(reader_function, inputs, reference_output):
    # test that the file readers work on a legitimate file,
    # and test that the output matches the reference output
    fn, offset = inputs
    actual_output = reader_function(fn, offset)
    assert_gem_results_equal(reference_output, _process_gemlog_data(actual_output))


# These tests ensure that it raises an exception when reading bad raw
# files. No need to classify the type of exception; any exception
# should be interpreted as a bad file.
@pytest.mark.parametrize('reader_function', [_read_with_cython,
                                             _read_with_pandas,
                                             _read_single])
def test__read_single_v0_9_empty(reader_function):
    with pytest.raises(EmptyRawFile):
        reader_function('../data/FILE0000.000')  # test empty file


@pytest.mark.parametrize('reader_function', [_read_with_cython,
                                             _read_with_pandas,
                                             _slow__read_single_v0_9])
def test__read_single_v0_9_corrupt(reader_function):
    with pytest.raises(Exception): # can't use CorruptRawFile: it's up to _read_single to make that call.
        reader_function('../data/FILE0023.096')  # test a malformed file

@pytest.mark.parametrize('reader_function', [_read_with_cython,
                                             _read_with_pandas,
                                             _slow__read_single_v0_9])
def test__read_single_v0_9_no_gps(reader_function):
    with pytest.raises(CorruptRawFileNoGPS):
        reader_function('../demo_missing_gps/raw_missing_gps/FILE0001.077')  # test a missing-gps file

###############################

## check that format v1.1 files can be read without error and that they give identical results to
## corresponding v0.91 files

def test_read_cython_v1_1_v0_91(): # at this point, only supporting cython reader
    ## initial file with a long GPS run at the beginning
    x = _read_with_cython('../data/v0.91/FILE0000.210')
    y = _read_with_cython('../data/v1.10/FILE0000.210')
    assert_gem_results_equal(_process_gemlog_data(x), _process_gemlog_data(y))

    ## ending file: almost full length, normal GPS runs
    x = _read_with_cython('../data/v0.91/FILE0001.210')
    y = _read_with_cython('../data/v1.10/FILE0001.210')
    assert_gem_results_equal(_process_gemlog_data(x), _process_gemlog_data(y))

    ## short initial file
    x = _read_with_cython('../data/v0.91/FILE0040.059')
    y = _read_with_cython('../data/v1.10/FILE0040.059')
    assert_gem_results_equal(_process_gemlog_data(x), _process_gemlog_data(y))


## incomplete last lines (power loss at the end of a recording) are ignored for AspenCSV0.2 and higher
def _make_aspen_fixture(src, dst, n_lines = 6000):
    """Write a smaller AspenCSV0.2 fixture so truncation tests run quickly."""
    with open(src, 'rb') as f_in, open(dst, 'wb') as f_out:
        for i, line in enumerate(f_in):
            if i >= n_lines:
                break
            f_out.write(line)


def _write_variant(src, dst, tail):
    """Copy src minus its last line (if incomplete) to dst, then append tail."""
    with open(src, 'rb') as f:
        lines = f.read().split(b'\n')
    # drop the last (incomplete or empty) element so only complete lines remain
    with open(dst, 'wb') as f:
        f.write(b'\n'.join(lines[:-1]) + b'\n' + tail)

def test_parser_drop_bad_last_line():
    _make_aspen_fixture('../data/AspenCSV/FILE0003.032', 'fixture.txt')
    _write_variant('fixture.txt', 'complete.txt', b'')
    ref = parse_gemfile(b'complete.txt', n_channels = 4, dt_ms = 5, drop_bad_last_line = True)
    n = len(ref[1])
    assert n > 0
    # complete last line is kept
    assert len(parse_gemfile(b'complete.txt', n_channels = 4, dt_ms = 5)[1]) == n
    tails = [b'D863,', b'D8', b'D', b'G,8171,70,2020,6,20,5,21,22.0,43.6', b'G,8171,70,2020,6,20,5,21,22.0,43.62226,-116.20594',
             b'H,773,3.77,3.784,-8.0,22.7', b'H,773,3.77,3.784,-8.0,22.7,31.9,3.29,6,1132,216,8375,101,0', b'#junk']
    for tail in tails:
        _write_variant('fixture.txt', 'bad.txt', tail) # no trailing newline
        x = parse_gemfile(b'bad.txt', n_channels = 4, dt_ms = 5, drop_bad_last_line = True)
        expected = n if tail == b'#junk' else n # unrecognized or incomplete trailing line never adds rows
        assert len(x[1]) == expected, tail
        assert np.array_equal(x[0], ref[0]) and np.array_equal(x[2], ref[2])
    # a truncated line that does end in a newline but can't be parsed is dropped too
    for tail in [b'G,8171,70,2020\n', b'H,773,3.77\n']:
        _write_variant('fixture.txt', 'bad.txt', tail)
        x = parse_gemfile(b'bad.txt', n_channels = 4, dt_ms = 5, drop_bad_last_line = True)
        assert len(x[1]) == n, tail
        # default behavior is unchanged: the bad line is kept
        assert len(parse_gemfile(b'bad.txt', n_channels = 4, dt_ms = 5)[1]) == n + 1
    # a bad line in the middle is not dropped
    with open('complete.txt', 'rb') as f:
        lines = f.read().split(b'\n')
    lines.insert(len(lines)//2, b'G,1,2')
    with open('mid.txt', 'wb') as f:
        f.write(b'\n'.join(lines))
    assert len(parse_gemfile(b'mid.txt', n_channels = 4, dt_ms = 5, drop_bad_last_line = True)[1]) == n + 1

def test_read_single_aspen_incomplete_last_line():
    _make_aspen_fixture('../data/AspenCSV/FILE0003.032', 'fixture.txt')
    _write_variant('fixture.txt', 'complete.txt', b'')
    ref = _read_single('complete.txt', version = 'AspenCSV0.2')
    for tail in [b'D863,', b'D8', b'G,8171,70,2020,6,20,5,21,22.0,43.6', b'H,773,3.7']:
        _write_variant('fixture.txt', 'bad.txt', tail)
        L = _read_single('bad.txt', version = 'AspenCSV0.2')
        assert np.array_equal(L['data'], ref['data'])
        assert len(L['gps']) == len(ref['gps'])
        assert len(L['metadata']) == len(ref['metadata'])

def test_aspen_drops_bad_last_line_versions():
    from gemlog.core import _aspen_drops_bad_last_line
    assert _aspen_drops_bad_last_line('AspenCSV0.2')
    assert _aspen_drops_bad_last_line('AspenCSV0.3')
    assert not _aspen_drops_bad_last_line('AspenCSV0.1')
    assert not _aspen_drops_bad_last_line('AspenCSV0.01')
    assert not _aspen_drops_bad_last_line('1.10')
    assert not _aspen_drops_bad_last_line(None)
