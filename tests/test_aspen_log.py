"""
Tests for reading log files from Aspen firmware 0.9 and later (AspenLog0.5), whose data are written directly as
miniSEED. The log text below follows the firmware's printLogHeader and printHealth (aspen-firmware/SD_functions.h and
aspen-firmware.ino).
"""
import os
import numpy as np
import pytest
import gemlog
from gemlog.core import _read_SN
from gemlog.aspen_log import read_aspen_log, extract_aspen_logs, find_aspen_logs, summarize_aspen_logs, clean_aspen_health

HEADER = """#AspenLog0.5
#Data are in the miniSEED file with the same name. This log holds everything else.
#T,tcxoSec,utcSec,sync: time of the next data sample (center of the ADC filter window); sync 0=no GPS time yet (utcSec blank), 1=GPS locked, 2=GPS holdover
#G,msPPS,msLag,yr,mo,day,hr,min,sec,lat,lon,z,HDOP,num_sats,age_ms
#H,tcxoSec,utcSec,batt,V,mA,tempC,RH,Vcc,maxSampleDiff_ms,maxReadTime_us,maxWriteTime_us,maxDrainTime_us,drainDiff_ms,gpsOn,sync,timeErr_us,tcxo_ppm,adc_ppm,maxRetimeTime_us,droppedRecords
#N,network,station,location,channel1,channel2,channel3,channel4 (SEED codes in the miniSEED file)
#miniSEED timing quality: 100 GPS locked, 90 GPS holdover, 80 recorded before the first GPS fix and timed from it, 0 no GPS time (time since startup, flagged questionable)
S,10024
#Firmware0.9.0-0_2026-10-18T00:00:00Z
C,128,64,64,64,5,1,15,20,0,0,1
N,XA,10024,,HDF,EHZ,EHN,EHE
"""


def health(tcxo, utc, sync):
    utc_text = '' if utc is None else f'{utc:.3f}'
    return f'H,{tcxo:.3f},{utc_text},3.96,2.568,-12.5,21.4,40.1,3.30,6,1154,59,1862,950,1,{sync},0.0,32.343,-97,324,0\n'


def gps(utc_s, lat = 43.62226):
    import datetime
    t = datetime.datetime.fromtimestamp(utc_s, tz = datetime.timezone.utc)
    return (f'G,1234.50,120.00,{t.year},{t.month},{t.day},{t.hour},{t.minute},{t.second:.3f},{lat:.6f},'
            '-116.205940,820.0,0.95,9,350\n')


UTC_LOCK = 1792201176.110  # 2026-10-17T01:39:36.110


def write_first_log(path, name = '20261017T013900_10024.log', incomplete_last_line = False):
    """A log that starts before the first GPS fix (blank UTC in early H lines) and then locks."""
    text = HEADER + 'T,0.134943,,0\n'
    text += 'R,137,$GNGGA,013856.211,,,,,0,0,,,M,,M,,*5D\n'
    for k in range(5):
        text += health(1.0 + k, None, 0)
    text += f'T,36.325378,{UTC_LOCK:.3f},1\n'
    for k in range(5):
        text += 'P,1234.50\n' + gps(int(UTC_LOCK) + 1 + k)
        text += health(37.0 + k, UTC_LOCK + 0.675 + k, 1)
    if incomplete_last_line:
        text += 'H,42.000,1792201182.0,3.9'
    (path / name).write_text(text)
    return path / name


def test_reads_health_gps_and_codes(tmp_path):
    L = read_aspen_log(str(write_first_log(tmp_path)))
    assert L['SN'] == '10024'
    assert L['codes'] == {'network': 'XA', 'station': '10024', 'location': '', 'channels': ['HDF', 'EHZ', 'EHN', 'EHE']}
    M = L['metadata']
    assert len(M) == 10
    assert list(M.columns[-1:]) == ['t']
    for col in ['batt', 'temp', 'V', 'mA', 'RH', 'tcxo_ppm', 'droppedRecords']:  # used by verify_huddle_test etc.
        assert col in M.columns
    # health lines before the first fix get UTC times from the T line at the fix
    np.testing.assert_allclose(M['t'].iloc[:5], UTC_LOCK - 36.325378 + np.arange(1.0, 6.0), atol = 1e-6)
    np.testing.assert_allclose(M['t'].iloc[5:], UTC_LOCK + 0.675 + np.arange(5), atol = 1e-6)
    G = L['gps']
    assert len(G) == 5
    np.testing.assert_allclose(G['t'], int(UTC_LOCK) + 1 + np.arange(5))
    np.testing.assert_allclose(G['lat'], 43.62226)


def test_reads_read_delays_when_present(tmp_path):
    # later firmware appends ADC2-4 read delays (us; blank for a channel that is off); older logs give NaN
    old = read_aspen_log(str(write_first_log(tmp_path)))['metadata']
    assert old[['readDelay2_us', 'readDelay3_us', 'readDelay4_us']].isna().all().all()
    text = HEADER.replace('droppedRecords\n', 'droppedRecords,readDelay2_us,readDelay3_us,readDelay4_us\n')
    text += f'T,36.325378,{UTC_LOCK:.3f},1\n' + health(37.0, UTC_LOCK + 0.675, 1).replace('\n', ',331,436,\n')
    (tmp_path / 'new.log').write_text(text)
    M = read_aspen_log(str(tmp_path / 'new.log'))['metadata']
    assert len(M) == 1
    assert M['droppedRecords'].iloc[0] == 0
    np.testing.assert_array_equal(M[['readDelay2_us', 'readDelay3_us']].iloc[0], [331, 436])
    assert np.isnan(M['readDelay4_us'].iloc[0])
    assert list(M.columns[-1:]) == ['t']


def test_ignores_incomplete_last_line(tmp_path):
    L = read_aspen_log(str(write_first_log(tmp_path, incomplete_last_line = True)))
    assert len(L['metadata']) == 10


def test_serial_number_found_like_other_raw_files(tmp_path):
    fn = write_first_log(tmp_path)
    assert _read_SN(str(fn)) == '10024'  # how gemconvert discovers serial numbers in a folder
    assert find_aspen_logs(str(tmp_path), '10024') == [str(fn)]
    assert find_aspen_logs(str(tmp_path), '10023') == []


def test_convert_writes_metadata_and_gps_tables(tmp_path):
    raw = tmp_path / 'raw'
    raw.mkdir()
    write_first_log(raw)
    second = HEADER + f'T,14436.325,{UTC_LOCK + 14400:.3f},2\n' + health(14437.0, UTC_LOCK + 14400.675, 2)
    second += gps(int(UTC_LOCK) + 14500, lat = 43.62227)
    (raw / '20261017T040000_10024.log').write_text(second)
    (raw / '20261017T013900_10024.mseed').write_bytes(b'\0' * 512)  # data files are left alone
    problems = gemlog.convert(str(raw), SN = '10024', convertedpath = str(tmp_path / 'converted'),
                              metadatapath = str(tmp_path / 'metadata'), gpspath = str(tmp_path / 'gps'))
    assert problems == []
    assert not os.path.exists(tmp_path / 'converted') or os.listdir(tmp_path / 'converted') == []
    import pandas as pd
    M = pd.read_csv(tmp_path / 'metadata' / '10024metadata_000.txt')
    assert len(M) == 11
    assert np.all(np.diff(M['t']) > 0)
    G = gemlog.read_gps(str(tmp_path / 'gps'), '10024')
    assert len(G) == 6
    assert list(G.columns) == ['year', 'date', 'lat', 'lon', 't']


def health_line(tcxo, utc, V = 4.02, mA = 9.0, temp = 20.0, RH = 50.0, batt = 4.0):
    return f'H,{tcxo:.3f},{utc:.3f},{batt},{V},{mA},{temp},{RH},3.30,6,1154,59,1862,950,0,1,0.0,32.343,-97,324,0\n'


def write_health_log(path, name, tcxo0, utc0, n, bad = {}):
    """n health lines, one per second starting at tcxo0; bad maps line index -> kwargs for health_line"""
    text = HEADER + f'T,{tcxo0:.3f},{utc0:.3f},1\n' + gps(int(utc0) + 1)
    for k in range(n):
        text += health_line(tcxo0 + k, utc0 + k, **bad.get(k, {}))
    (path / name).write_text(text)


def test_clean_health_removes_bad_values(tmp_path):
    bad = {1: dict(temp = 130.0, RH = 3.3), 2: dict(V = -32.76, mA = 409.6), 3: dict(temp = -44.8), 4: dict(mA = -2662.4)}
    write_health_log(tmp_path, 'b_10024.log', 5000, UTC_LOCK, 6, bad)
    M = read_aspen_log(str(tmp_path / 'b_10024.log'))['metadata']
    C, removed = clean_aspen_health(M)
    assert C['temp'].isna().tolist() == [False, True, False, True, False, False]
    assert C['RH'].isna().tolist() == [False, True, False, False, False, False]
    assert C['V'].isna().tolist() == [False, False, True, False, False, False]
    assert C['mA'].isna().tolist() == [False, False, True, False, True, False]
    assert C['batt'].notna().all()
    # first 960 s after startup: power monitor, temperature, and RH are all removed, battery (analog) is kept
    early = pd_read(tmp_path, 'c_10024.log', 100)
    C, _ = clean_aspen_health(early)
    assert C[['V', 'mA', 'temp', 'RH']].isna().all().all()
    assert C['batt'].notna().all()


def pd_read(path, name, tcxo0):
    write_health_log(path, name, tcxo0, UTC_LOCK, 4)
    return read_aspen_log(str(path / name))['metadata']


def test_summary_figures_and_no_warnings(tmp_path):
    write_health_log(tmp_path, '20261017T000000_10024.log', 5000, UTC_LOCK, 100,
                     {50: dict(temp = 130.0, RH = 3.3)})
    write_health_log(tmp_path, '20261017T000200_10024.log', 5100, UTC_LOCK + 100, 100, {0: dict(temp = 30.0)})
    R = summarize_aspen_logs(str(tmp_path), verbose = False)
    S = R['summaries']['10024']
    assert len(S['files']) == 2
    assert S['duration_hours'] == pytest.approx(199 / 3600)
    assert S['stats'].loc['temp', 'max'] == 30.0
    assert S['stats'].loc['temp', 'n'] == 199
    assert S['stats'].loc['V', 'median'] == pytest.approx(4.02)
    assert S['stats'].loc['RH', 'pct_removed'] == pytest.approx(0.5)
    assert S['warnings'] == []


def test_summary_warnings(tmp_path):
    write_health_log(tmp_path, '20261017T000000_10024.log', 5000, UTC_LOCK, 100, {5: dict(batt = 3.2)})
    write_health_log(tmp_path, '20261017T001000_10024.log', 100, UTC_LOCK + 600, 100)  # restart, 500 s gap
    with open(tmp_path / '20261017T001000_10024.log', 'a') as f:
        f.write('H,300.000,17922')  # incomplete last line
    (tmp_path / 'notes.log').write_text('not an Aspen log\n')
    R = summarize_aspen_logs(str(tmp_path), verbose = False)
    W = ' | '.join(R['summaries']['10024']['warnings'])
    assert 'battery voltage dropped to 3.20 V' in W
    assert 'logger restarted' in W
    assert 'missing health data' in W
    assert 'incomplete line' in W
    assert any('notes.log' in p for p in R['problems'])


def test_summary_figure_is_saved(tmp_path):
    import matplotlib
    matplotlib.use('Agg')
    write_health_log(tmp_path, '20261017T000000_10024.log', 5000, UTC_LOCK, 100, {50: dict(temp = 130.0, RH = 3.3)})
    out = tmp_path / 'fig.png'
    R = summarize_aspen_logs(str(tmp_path), verbose = False, plot = str(out))
    assert out.stat().st_size > 1000
    fig = gemlog.plot_aspen_summary(R)
    assert len(fig.axes) == 4
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_clean_health_dead_sensors_and_high_voltage_supply(tmp_path):
    bad = {1: dict(temp = -45.0, RH = 0.0), 2: dict(V = -0.01, mA = -0.8), 3: dict(V = 3.2, mA = 7.0)}
    write_health_log(tmp_path, 'a_10024.log', 5000, UTC_LOCK, 5, bad)
    C, _ = clean_aspen_health(read_aspen_log(str(tmp_path / 'a_10024.log'))['metadata'])
    assert C['temp'].isna().tolist() == [False, True, False, False, False]
    assert C['RH'].isna().tolist() == [False, True, False, False, False]
    assert C['V'].isna().tolist() == [False, False, True, True, False]  # row 3: V does not match batt (4.0)
    assert C['mA'].isna().tolist() == [False, False, True, True, False]
    # a logger running from a 5.4 V supply is not bad data
    write_health_log(tmp_path, 'b_10024.log', 5000, UTC_LOCK, 5, {k: dict(V = 5.4, batt = 5.45) for k in range(5)})
    C, _ = clean_aspen_health(read_aspen_log(str(tmp_path / 'b_10024.log'))['metadata'])
    assert C[['batt', 'V', 'mA']].notna().all().all()


def test_plot_true_shows_figures(tmp_path, monkeypatch):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.close('all')
    shown = []
    monkeypatch.setattr(plt, 'show', lambda *a, **k: shown.append(len(plt.get_fignums())))
    write_health_log(tmp_path, '20261017T000000_10024.log', 5000, UTC_LOCK, 100)
    summarize_aspen_logs(str(tmp_path), verbose = False, plot = True)
    assert shown == [1]
    plt.close('all')
