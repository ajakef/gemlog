"""
Tests for AspenCSV0.4, where the logger writes samples already evenly spaced on the UTC grid and T lines give their
times. The raw files are synthetic (see aspen04_files.py).
"""
import os
import numpy as np
import obspy
import pytest
import gemlog
from gemlog.core import read_gem, _read_aspen04_single, _read_format_version, _read_config_aspen, _read_SN
from gemlog.exceptions import CorruptRawFileInadequateGPS
from gemlog.parsers import parse_aspen04
from aspen04_files import Aspen04Writer, health_line, gps_line

DT_MS = 5
UTC0_MS = 1792108800000  # 2026-10-16 00:00:00 UTC
BOOT_UTC_MS = UTC0_MS + 12347  # true UTC when the logger's TCXO clock read 0 (whole ms keeps the test exact;
                               # not a multiple of 5 ms, so samples before the first fix are off the grid)


def signal(utc_s, ch):
    """Test signal in counts: a slow wave (mostly compressed lines) plus a 40 Hz part (some uncompressed lines)."""
    t = np.asarray(utc_s) - UTC0_MS / 1000
    return np.round(100000 * ch + 3000 * np.sin(2 * np.pi * 3 * t) + 300 * np.sin(2 * np.pi * 40 * t + ch))


def make_recording(directory, n_files = 2, samples_per_file = 6000, unsynced = 1000, gap_at = None,
                   gps = True, SN = '00022'):
    """
    Writes raw files like a logger that starts without GPS time, gets its first fix after `unsynced` samples, and
    then stays locked. The TCXO is perfect here, so the logger's own times are exact.
    Returns a list of (utc_ms, values) for every sample written, where utc_ms is the true sample time.
    """
    written = []
    tcxo_ms = 200  # first sample
    locked = False
    utc_ms = None
    writer = None
    file_num = 0
    count_in_file = samples_per_file  # forces a new file at the start
    k = 0
    total = n_files * samples_per_file
    pending_anchor = False
    while k < total:
        if count_in_file == samples_per_file:
            if writer is not None:
                writer.close()
            writer = Aspen04Writer(os.path.join(directory, f'FILE{file_num:04d}.{SN[-3:]}'), SN = SN)
            file_num += 1
            count_in_file = 0
            pending_anchor = True
        if (not locked) and gps and (k == unsynced):
            # first GPS fix: jump to the first UTC grid point at or after the next nominal sample time
            nominal = BOOT_UTC_MS + tcxo_ms
            utc_ms = -(-nominal // DT_MS) * DT_MS
            tcxo_ms = utc_ms - BOOT_UTC_MS
            locked = True
            pending_anchor = True
        if (gap_at is not None) and (k == gap_at):
            tcxo_ms += 16 * DT_MS  # the logger missed some readings
            utc_ms += 16 * DT_MS
            pending_anchor = True
        true_utc_ms = BOOT_UTC_MS + tcxo_ms
        if pending_anchor:
            writer.anchor(tcxo_ms / 1000, utc_ms / 1000 if locked else None, 1 if locked else 0)
            pending_anchor = False
        values = [signal(true_utc_ms / 1000, ch) for ch in range(4)]
        writer.data(values)
        written.append((true_utc_ms, values))
        if k % 200 == 35:
            writer.line(health_line((tcxo_ms * 1024 // 1000) % 8192, 1 if locked else 0))
        if locked and (k % 200 == 100):
            pps_s = true_utc_ms // 1000
            writer.line(f'P,{((pps_s * 1000 - BOOT_UTC_MS) * 1.024) % 8192:.2f}')
            writer.line(gps_line(pps_s, ((pps_s * 1000 - BOOT_UTC_MS) * 1.024) % 8192))
        tcxo_ms += DT_MS
        if locked:
            utc_ms += DT_MS
        k += 1
        count_in_file += 1
    writer.close()
    return written


def test_header_is_read(tmp_path):
    make_recording(str(tmp_path), n_files = 1, samples_per_file = 1500)
    fn = str(tmp_path / 'FILE0000.022')
    assert _read_format_version(fn) == 'AspenCSV0.4'
    assert _read_SN(fn) == '00022'
    config = _read_config_aspen(fn)
    assert config['sample_int_ms'] == 5 and config['n_channels'] == 4


def test_parser_keeps_raw_values_and_drops_incomplete_last_line(tmp_path):
    fn = tmp_path / 'FILE0000.022'
    fn.write_text('#AspenCSV0.4\nS,00022\nC,128,32,0,0,5,1,15,25,0,0,1\n'
                  'T,1.000000,,0\n1000,-2000\nnl\nm\nH,1,2,3\n12,-30\nT,2.000000,1792108800.005,1\n'
                  '5,6\nG,1,2\nnm\n7')  # the last line has no newline (power lost while writing)
    values, others = parse_aspen04(str(fn).encode(), n_channels = 2, n_row = 100)
    np.testing.assert_array_equal(values, [[1000, -2000], [1, -1], [0, 0], [12, -30], [5, 6], [1, 0]])
    assert [(i, line.split(b',')[0]) for i, line in others] == [(0, b'T'), (3, b'H'), (4, b'T'), (5, b'G')]


def test_values_and_times_rebuilt_from_t_lines(tmp_path):
    fn = tmp_path / 'FILE0000.022'
    fn.write_text('#AspenCSV0.4\nS,00022\nC,128,32,0,0,5,1,15,25,0,0,1\n'
                  'T,1.000000,,0\n1000,-2000\nnl\nm\n'
                  'T,1.020000,1792108800.005,1\n5,6\nnm\n')
    F = _read_aspen04_single(str(fn))
    np.testing.assert_array_equal(F['values'], [[1000, -2000], [1001, -2001], [1001, -2001], [5, 6], [6, 6]])
    np.testing.assert_array_equal(F['sync'], [0, 0, 0, 1, 1])
    assert np.all(np.isnan(F['utc'][:3]))
    np.testing.assert_allclose(F['utc'][3:], [1792108800.005, 1792108800.010], rtol = 0, atol = 1e-6)
    np.testing.assert_allclose(F['tcxo'], [1.0, 1.005, 1.010, 1.020, 1.025], rtol = 0, atol = 1e-9)


def test_read_gem_gives_exact_samples_on_utc_grid(tmp_path):
    written = make_recording(str(tmp_path), n_files = 2, samples_per_file = 6000, unsynced = 1000)
    L = read_gem(str(tmp_path), SN = '00022')
    st = L['data']
    assert sorted(set(tr.stats.channel for tr in st)) == ['HDF', 'XX1', 'XX2', 'XX3']
    tr_hdf = st.select(channel = 'HDF')
    # before the first fix (interpolated onto the grid) and after it: two traces, despite the file boundary
    assert len(tr_hdf) == 2
    locked = tr_hdf[1]
    assert locked.stats.delta == 0.005
    assert locked.stats.starttime.ns % 5000000 == 0  # on the 5 ms grid
    locked_written = written[1000:]
    assert locked.stats.npts == len(locked_written)
    assert locked.stats.starttime.ns == locked_written[0][0] * 1000000
    np.testing.assert_array_equal(locked.data, [v[0] for _, v in locked_written])
    xx3 = st.select(channel = 'XX3')[1]
    np.testing.assert_array_equal(xx3.data, [v[3] for _, v in locked_written])

    # samples recorded before the first fix get their time from the fix and are interpolated onto the grid
    early = tr_hdf[0]
    assert early.stats.starttime.ns % 5000000 == 0
    t_early = early.stats.starttime.timestamp + np.arange(early.stats.npts) * 0.005
    assert t_early[0] >= written[0][0] / 1000 and t_early[-1] <= written[999][0] / 1000
    assert np.max(np.abs(early.data - signal(t_early, 0))) < 20  # cubic spline error for the 40 Hz part

    # metadata and gps
    assert np.all(np.diff(L['metadata']['t']) > 0)
    assert L['metadata']['t'].iloc[0] >= written[0][0] / 1000 - 1
    assert 'sync' in L['metadata'] and 'adc_ppm' in L['metadata']
    assert len(L['gps']) > 0


def test_gap_in_samples_splits_trace(tmp_path):
    written = make_recording(str(tmp_path), n_files = 1, samples_per_file = 6000, unsynced = 0, gap_at = 3000)
    st = read_gem(str(tmp_path), SN = '00022')['data'].select(channel = 'HDF')
    assert len(st) == 2
    assert st[0].stats.npts == 3000 and st[1].stats.npts == 3000
    assert (st[1].stats.starttime - st[0].stats.endtime) == pytest.approx(17 * 0.005)
    np.testing.assert_array_equal(st[1].data, [v[0] for _, v in written[3000:]])


def test_no_gps_time(tmp_path):
    make_recording(str(tmp_path), n_files = 1, samples_per_file = 2000, gps = False)
    with pytest.raises(CorruptRawFileInadequateGPS):
        read_gem(str(tmp_path), SN = '00022')
    # without requiring GPS, times count from when the logger started (as in other formats without GPS)
    st = read_gem(str(tmp_path), SN = '00022', require_gps = False)['data'].select(channel = 'HDF')
    assert st[0].stats.npts >= 1990
    assert st[0].stats.starttime < obspy.UTCDateTime(1000)


def test_convert_writes_exact_data(tmp_path):
    raw = tmp_path / 'raw'
    raw.mkdir()
    written = make_recording(str(raw), n_files = 2, samples_per_file = 6000, unsynced = 0)
    gemlog.convert(str(raw), convertedpath = str(tmp_path / 'converted'), metadatapath = str(tmp_path / 'metadata'),
                   gpspath = str(tmp_path / 'gps'), SN = '00022', units = 'counts')
    st = obspy.read(str(tmp_path / 'converted' / '*'))
    tr = st.select(channel = 'HDF').merge()[0]
    first = (tr.stats.starttime.ns // 1000000 - written[0][0]) // DT_MS  # convert starts on a whole second
    assert 0 <= first < 200
    assert tr.stats.starttime.ns == written[first][0] * 1000000
    np.testing.assert_array_equal(tr.data, [v[0] for _, v in written[first:first + tr.stats.npts]])
    assert first + tr.stats.npts >= len(written) - 1
