"""
Log files from Aspen firmware 0.9 and later (format AspenLog0.5).

These loggers write their data directly as miniSEED (in counts, already on the UTC sample grid), so the data need no
conversion. Everything else goes in a text log with the same name: health (H) lines, GPS (G) lines, and time (T)
lines. This module turns the logs into the same metadata and GPS tables that gemlog.convert writes for older formats,
so tools like verify_huddle_test, read_gps, and summarize_gps work the same way.
"""
import glob
import os
import numpy as np
import pandas as pd
from gemlog.core import _gps_in_bounds, _make_gps_time, _reformat_GPS, _make_filename
from gemlog.exceptions import CorruptRawFile

_LOG_PREFIX = '#AspenLog'
_H_COLS = ['tcxoSec', 'utcSec', 'batt', 'V', 'mA', 'temp', 'RH', 'Vcc', 'maxSampleDiff', 'maxReadTime',
           'maxWriteTime', 'maxDrainTime', 'drainDiff', 'gpsOnFlag', 'sync', 'timeErr_us', 'tcxo_ppm', 'adc_ppm',
           'maxRetimeTime', 'droppedRecords']
## H columns added by later firmware (after 0.8.6-1): NaN for logs without them. readDelay: us from ADC1's DRDY to
## reading ADC2-4 (blank if the channel is off), which sets the firmware's correction for those readings' age.
_H_OPTIONAL_COLS = ['readDelay2_us', 'readDelay3_us', 'readDelay4_us']
## Health (H) line values outside these ranges are not physically possible for the sensors, so they are bad reads
_VALID_RANGE = {'batt': (2.5, 16), 'V': (2.5, 16), 'mA': (-200, 200), 'temp': (-40, 85), 'RH': (0, 100)}
_SUMMARY_COLS = ['batt', 'V', 'mA', 'temp', 'RH']
_I2C_COLS = ['V', 'mA', 'temp', 'RH']  # read over I2C, so affected by the startup garbage (not batt, an analog read)
_MAX_BATT_V_DISAGREE = 0.5  # V (power monitor) must match batt (analog) within this, or V and mA are bad
_BOOT_GARBAGE_SEC = 960  # firmware 0.7.0-0.8.2 give garbage I2C values for ~900-910 s after every startup
## limits for warnings in summarize_aspen_logs
_MAX_H_GAP_SEC = 10  # health lines are written every second
_MAX_FILE_GAP_SEC = 30
_MAX_TIME_ERR_US = 5000
_MAX_BAD_FRACTION = 0.01  # of values after startup
_MAX_NO_GPS_FRACTION = 0.1
_MAX_BATT_DISAGREE_V = 0.3
_LOW_BATT_V = 3.4
_LOW_VCC = 3.2
_TEMP_OK = (-30, 60)
_HIGH_RH = 95
_G_COLS = ['msPPS', 'msLag', 'year', 'month', 'day', 'hour', 'minute', 'second', 'lat', 'lon']


def _is_aspen_log(filename):
    try:
        with open(filename, errors = 'ignore') as f:
            return f.readline().startswith(_LOG_PREFIX)
    except OSError:
        return False


def read_aspen_log(filename):
    """
    Read one Aspen log file (AspenLog0.5).

    Parameters
    ----------
    filename : str
        Path of the log file.

    Returns
    -------
    dict with keys:
        - SN : str, serial number
        - codes : dict of SEED codes used in the miniSEED file (network, station, location, channels)
        - metadata : pandas.DataFrame, one row per health (H) line; 't' is UTC (seconds since 1970). Columns added by
          later firmware (_H_OPTIONAL_COLS) are always present, NaN when the log doesn't have them.
        - gps : pandas.DataFrame, one row per valid GPS fix (G line); 't' is UTC (seconds since 1970)
        - firmware : str, firmware version line ('' if missing)
        - n_corrupt_lines : int, number of H, G, T, S, or N lines that could not be read
        - truncated : bool, True if the file ends in an incomplete line (power lost while writing)
    """
    if not _is_aspen_log(filename):
        raise CorruptRawFile(f'{filename} is not an Aspen log file')
    SN = ''
    codes = {}
    firmware = ''
    n_corrupt = 0
    truncated = False
    T_rows, H_rows, G_rows = [], [], []
    with open(filename, errors = 'ignore') as f:
        for line in f:
            if not line.endswith('\n'):
                truncated = True
                break  # incomplete last line (power lost while writing)
            if line.startswith('#Firmware'):
                firmware = line.strip()[len('#Firmware'):]
            fields = line.rstrip('\r\n').split(',')
            try:
                if fields[0] == 'T':
                    T_rows.append([float(fields[1]), float(fields[2]) if fields[2] != '' else np.nan])
                elif fields[0] == 'H':
                    n = len(_H_COLS) + len(_H_OPTIONAL_COLS)
                    vals = [float(x) if x != '' else np.nan for x in fields[1:n+1]]
                    if len(vals) >= len(_H_COLS):
                        H_rows.append(vals + [np.nan] * (n - len(vals)))
                    else:
                        n_corrupt += 1
                elif fields[0] == 'G':
                    vals = [float(x) for x in fields[1:len(_G_COLS)+1]]
                    if len(vals) == len(_G_COLS):
                        G_rows.append(vals)
                    else:
                        n_corrupt += 1
                elif fields[0] == 'S':
                    SN = fields[1]
                elif fields[0] == 'N':
                    codes = {'network': fields[1], 'station': fields[2], 'location': fields[3],
                             'channels': fields[4:8]}
            except (ValueError, IndexError):
                n_corrupt += 1  # corrupt line: skip it

    ## health lines: UTC time, or (before the logger's first GPS fix) time since startup, converted using the first
    ## T line in this log that has both
    M = pd.DataFrame(H_rows, columns = _H_COLS + _H_OPTIONAL_COLS)
    T = np.array(T_rows, dtype = float).reshape(-1, 2)
    timed = ~np.isnan(T[:,1])
    if np.any(timed) and M['utcSec'].isna().any():
        offset = T[timed][0, 1] - T[timed][0, 0]
        untimed = M['utcSec'].isna()
        M.loc[untimed, 'utcSec'] = M.loc[untimed, 'tcxoSec'] + offset
    M = M.rename(columns = {'utcSec': 't'})
    M = M[[c for c in M.columns if c != 't'] + ['t']]  # 't' last, as in metadata from other formats

    if len(G_rows) > 0:
        G = pd.DataFrame(G_rows, columns = _G_COLS)
        G = G.loc[_gps_in_bounds(G), :].copy()
        G['t'] = [float(_make_gps_time(row)) for row in G[_G_COLS].values]
        G = G.loc[~G['t'].isna(), :].reset_index(drop = True)
    else:
        G = pd.DataFrame(columns = _G_COLS + ['t'])
    return {'SN': SN, 'codes': codes, 'metadata': M, 'gps': G, 'firmware': firmware,
            'n_corrupt_lines': n_corrupt, 'truncated': truncated}


def find_aspen_logs(rawpath, SN = None):
    """List the Aspen log files in rawpath (for serial number SN, if given), sorted by name (start time)."""
    logs = []
    for fn in sorted(glob.glob(os.path.join(rawpath, '*.log'))):
        if not _is_aspen_log(fn):
            continue
        if SN is not None:
            try:
                if read_aspen_log_SN(fn) != str(SN):
                    continue
            except CorruptRawFile:
                continue
        logs.append(fn)
    return logs


def read_aspen_log_SN(filename):
    """Serial number from an Aspen log's S line."""
    with open(filename, errors = 'ignore') as f:
        for i, line in enumerate(f):
            if line.startswith('S,'):
                return line.rstrip('\r\n').split(',')[1]
            if i > 50:
                break
    raise CorruptRawFile(f'{filename}: missing serial number')


def extract_aspen_logs(rawpath = '.', SN = None, metadatapath = 'metadata', gpspath = 'gps'):
    """
    Write metadata (health) and GPS tables from Aspen log files, in the same form gemlog.convert writes them for
    older formats: metadatapath/<SN>metadata_NNN.txt and gpspath/<SN>gps_NNN.txt.

    Parameters
    ----------
    rawpath : str
        Folder containing the logger's .log (and .mseed) files.
    SN : str, optional
        Serial number to process. By default, all serial numbers found in rawpath.
    metadatapath, gpspath : str
        Output folders (created if needed).

    Returns
    -------
    list of str: problems found (e.g. unreadable files)
    """
    problems = []
    logs = find_aspen_logs(rawpath, SN)
    by_SN = {}
    for fn in logs:
        try:
            L = read_aspen_log(fn)
        except CorruptRawFile as e:
            problems.append(str(e))
            continue
        by_SN.setdefault(L['SN'], []).append(L)
    for sn, Ls in by_SN.items():
        print(f'{sn}: {len(Ls)} Aspen log files; the data are already miniSEED (counts), so only metadata and GPS '
              'are extracted')
        for path in [metadatapath, gpspath]:
            os.makedirs(path, exist_ok = True)
        M = pd.concat([L['metadata'] for L in Ls], ignore_index = True).sort_values('t', kind = 'stable')
        M.to_csv(_make_filename(metadatapath, sn, 'metadata'), index = False)
        gps_chunks = [L['gps'] for L in Ls if L['gps'].shape[0] > 0]
        if len(gps_chunks) > 0:
            G = pd.concat(gps_chunks, ignore_index = True).sort_values('t', kind = 'stable')
            _reformat_GPS(G).to_csv(_make_filename(gpspath, sn, 'gps'), index = False)
        else:
            problems.append(f'{sn}: no GPS fixes in the log files')
    return problems


def clean_aspen_health(M):
    """
    Set obviously bad values in the health table to NaN.

    Bad values are: batt, V, mA, temp, and RH outside physically possible ranges; temp = 130, temp = -45, RH = 3.3, or
    RH = 0 (the sensor read failed, and both values are bad); V and mA when V is out of range or differs from batt by
    more than 0.5 V (same chip); and V, mA, temp, and RH in the first 960 s after startup (known
    power-monitor/temperature/RH garbage in firmware 0.7.0-0.8.2).

    Parameters
    ----------
    M : pandas.DataFrame
        Health table from read_aspen_log (needs tcxoSec and the columns batt, V, mA, temp, RH).

    Returns
    -------
    (pandas.DataFrame, pandas.DataFrame): copy of M with bad values set to NaN, and a boolean table (same index)
    saying which values were removed, with columns batt, V, mA, temp, RH.
    """
    C = M.copy()
    bad = pd.DataFrame(False, index = M.index, columns = _SUMMARY_COLS)
    for col, (lo, hi) in _VALID_RANGE.items():
        bad[col] |= ~M[col].between(lo, hi) & M[col].notna()
    ## temp/RH sensor: 130 and 3.3 are failed reads; -45 and 0 are what it returns when it does not respond at all
    sensor_failed = (M['temp'] == 130) | (M['RH'] == 3.3) | (M['temp'] == -45) | (M['RH'] == 0)
    bad['temp'] |= sensor_failed
    bad['RH'] |= sensor_failed
    ## V and mA come from the same chip: if V is impossible or does not match the battery measured separately (batt),
    ## the chip's readings are bad
    power_failed = bad['V'] | (bad['batt'].eq(False) & ((M['V'] - M['batt']).abs() > _MAX_BATT_V_DISAGREE))
    bad['V'] |= power_failed
    bad['mA'] |= power_failed
    for col in _I2C_COLS:
        bad[col] |= (M['tcxoSec'] < _BOOT_GARBAGE_SEC)
    for col in _SUMMARY_COLS:
        C.loc[bad[col], col] = np.nan
    return C, bad


def _utc_text(t):
    return 'unknown' if not np.isfinite(t) else pd.to_datetime(t, unit = 's').strftime('%Y-%m-%d %H:%M:%S')


def summarize_aspen_logs(path = '.', SN = None, verbose = True, plot = None):
    """
    Read several Aspen log files and summarize them, with warnings about unusual or problematic behavior.

    Parameters
    ----------
    path : str or list of str
        Folder containing the .log files, or a list of log file names.
    SN : str, optional
        Serial number to summarize. By default, all serial numbers found.
    verbose : bool
        If True, print the report.
    plot : bool or str, optional
        True shows a summary figure for each serial number (see plot_aspen_summary). A string saves the figures
        instead: it is a file name or a folder; a folder gets <SN>_health_summary.png.

    Returns
    -------
    dict with keys:
        - summaries : dict, one entry per serial number, each a dict with
            - files : list of file names used (in time order)
            - start, end : UTC text of the first and last health line
            - duration_hours : float
            - firmware : list of firmware versions found
            - codes : SEED codes (from the first file)
            - n_gps_fixes : int
            - stats : pandas.DataFrame with one row each for batt, V, mA, temp, RH and columns n, min, median, mean,
              max, pct_removed. Obviously bad values are left out (see clean_aspen_health).
            - warnings : list of str
            - health : pandas.DataFrame, health lines with bad values set to NaN (used for the figure)
        - problems : list of str, files that could not be used
    """
    problems = []
    if isinstance(path, (str, os.PathLike)):
        names = sorted(glob.glob(os.path.join(path, '*.log')))
    else:
        names = list(path)
    logs = []
    for fn in names:
        try:
            L = read_aspen_log(fn)
        except (CorruptRawFile, OSError):
            problems.append(f'{fn}: skipped, not a readable Aspen log file')
            continue
        L['filename'] = fn
        if SN is None or L['SN'] == str(SN):
            logs.append(L)

    by_SN = {}
    for L in logs:
        by_SN.setdefault(L['SN'], []).append(L)
    summaries = {}
    for sn, Ls in by_SN.items():
        summaries[sn] = _summarize_one_logger(sn, Ls)
    result = {'summaries': summaries, 'problems': problems}
    if verbose:
        print(format_aspen_summary(result))
    if plot is True:
        import matplotlib.pyplot as plt
        for sn in summaries:
            plot_aspen_summary(result, sn)
        plt.show()
    elif plot:
        import matplotlib.pyplot as plt
        for sn in summaries:
            fn = os.path.join(plot, f'{sn}_health_summary.png') if os.path.isdir(plot) else plot
            if len(summaries) > 1 and not os.path.isdir(plot):
                root, ext = os.path.splitext(plot)
                fn = f'{root}_{sn}{ext}'
            plt.close(plot_aspen_summary(result, sn, filename = fn))
    return result


def _summarize_one_logger(sn, Ls):
    warnings = []
    ## time order: by first UTC time if known, otherwise by file name
    first_t = [L['metadata']['t'].min() if len(L['metadata']) > 0 else np.nan for L in Ls]
    order = sorted(range(len(Ls)), key = lambda i: (np.isnan(first_t[i]), first_t[i], Ls[i]['filename']))
    Ls = [Ls[i] for i in order]

    ## per-file checks
    for L in Ls:
        name = os.path.basename(L['filename'])
        M = L['metadata']
        if len(M) == 0:
            warnings.append(f'{name}: no health lines (empty log)')
            continue
        if L['truncated']:
            warnings.append(f'{name}: ends in an incomplete line (power lost or card removed while writing)')
        if L['n_corrupt_lines'] > 0:
            warnings.append(f"{name}: {L['n_corrupt_lines']} unreadable lines skipped")
        dt = np.diff(M['tcxoSec'].values)
        if np.any(dt > _MAX_H_GAP_SEC):
            warnings.append(f'{name}: {np.sum(dt > _MAX_H_GAP_SEC)} gaps of more than {_MAX_H_GAP_SEC} s between '
                            f'health lines (longest {dt.max():.0f} s); logger may have stalled')
        if np.any(dt < 0):
            warnings.append(f'{name}: time since startup runs backward in {np.sum(dt < 0)} places')
        if M['droppedRecords'].max() > 0:
            warnings.append(f"{name}: {M['droppedRecords'].max():.0f} data records dropped (droppedRecords)")
        err = M.loc[M['sync'] >= 1, 'timeErr_us'].abs().max()
        if err > _MAX_TIME_ERR_US:
            warnings.append(f'{name}: clock error reached {err:.0f} us while GPS time was known '
                            f'(limit {_MAX_TIME_ERR_US} us)')
        if (M['sync'] == 0).all():
            warnings.append(f'{name}: no GPS time during the whole file')
        elif L['gps'].shape[0] == 0:
            warnings.append(f'{name}: has GPS time but no valid GPS fixes')

    ## checks between files
    usable = [L for L in Ls if len(L['metadata']) > 0]
    for a, b in zip(usable[:-1], usable[1:]):
        na, nb = os.path.basename(a['filename']), os.path.basename(b['filename'])
        Ma, Mb = a['metadata'], b['metadata']
        if Mb['tcxoSec'].iloc[0] < Ma['tcxoSec'].iloc[-1]:
            warnings.append(f'{nb}: logger restarted (time since startup is lower than at the end of {na})')
        gap = Mb['t'].iloc[0] - Ma['t'].iloc[-1]
        if np.isfinite(gap) and gap > _MAX_FILE_GAP_SEC:
            warnings.append(f'{na} -> {nb}: {gap:.0f} s of missing health data between the files')
        elif np.isfinite(gap) and gap < -1:
            warnings.append(f'{na} -> {nb}: files overlap in time by {-gap:.0f} s')
    if len({L['firmware'] for L in Ls}) > 1:
        warnings.append('more than one firmware version: ' + ', '.join(sorted({L['firmware'] for L in Ls})))
    if len({str(L['codes']) for L in Ls}) > 1:
        warnings.append('SEED codes (N line) differ between files')

    ## summary figures, from all files together
    M = pd.concat([L['metadata'] for L in usable], ignore_index = True) if usable else pd.DataFrame()
    stats = pd.DataFrame(index = _SUMMARY_COLS, columns = ['n', 'min', 'median', 'mean', 'max', 'pct_removed'],
                         dtype = float)
    t = M['t'] if len(M) > 0 else pd.Series(dtype = float)
    if len(M) > 0:
        C, bad = clean_aspen_health(M)
        for col in _SUMMARY_COLS:
            x = C[col].dropna()
            stats.loc[col] = [len(x), x.min(), x.median(), x.mean(), x.max(), 100 * bad[col].mean()]
        after_boot = M['tcxoSec'] >= _BOOT_GARBAGE_SEC
        for col in _SUMMARY_COLS:
            if after_boot.sum() > 0 and bad.loc[after_boot, col].mean() > _MAX_BAD_FRACTION:
                warnings.append(f'{col}: {100 * bad.loc[after_boot, col].mean():.1f}% of values after the first '
                                f'{_BOOT_GARBAGE_SEC} s were bad and removed from the summary')
        ## limits on the cleaned values
        if stats.loc['batt', 'n'] > 0 and stats.loc['batt', 'min'] < _LOW_BATT_V:
            warnings.append(f"battery voltage dropped to {stats.loc['batt', 'min']:.2f} V (below {_LOW_BATT_V} V)")
        if stats.loc['V', 'n'] > 0 and stats.loc['V', 'min'] < _LOW_BATT_V:
            warnings.append(f"power-monitor voltage dropped to {stats.loc['V', 'min']:.2f} V (below {_LOW_BATT_V} V)")
        if stats.loc['temp', 'n'] > 0 and not _TEMP_OK[0] <= stats.loc['temp', 'min'] <= stats.loc['temp', 'max'] <= _TEMP_OK[1]:
            warnings.append(f"temperature range {stats.loc['temp', 'min']:.1f} to {stats.loc['temp', 'max']:.1f} C "
                            f'is outside {_TEMP_OK[0]} to {_TEMP_OK[1]} C')
        if stats.loc['RH', 'n'] > 0 and stats.loc['RH', 'max'] > _HIGH_RH:
            warnings.append(f"relative humidity reached {stats.loc['RH', 'max']:.0f}% (condensation risk)")
        vcc = M['Vcc'].dropna()
        if len(vcc) > 0 and vcc.min() < _LOW_VCC:
            warnings.append(f'processor voltage (Vcc) dropped to {vcc.min():.2f} V (below {_LOW_VCC} V)')
        if (M['sync'] == 0).mean() > _MAX_NO_GPS_FRACTION:
            warnings.append(f"no GPS time for {100 * (M['sync'] == 0).mean():.0f}% of the health lines")
    n_gps = int(sum(L['gps'].shape[0] for L in Ls))
    if n_gps == 0:
        warnings.append('no valid GPS fixes in any file')
    return {'files': [L['filename'] for L in Ls],
            'start': _utc_text(t.min()) if len(t) else 'unknown',
            'end': _utc_text(t.max()) if len(t) else 'unknown',
            'duration_hours': float((t.max() - t.min()) / 3600) if len(t) and t.notna().any() else np.nan,
            'firmware': sorted({L['firmware'] for L in Ls if L['firmware']}),
            'codes': Ls[0]['codes'],
            'n_gps_fixes': n_gps,
            'stats': stats,
            'warnings': warnings,
            'health': C if len(M) > 0 else pd.DataFrame()}


def format_aspen_summary(result):
    """Text report from the output of summarize_aspen_logs."""
    lines = []
    for sn, S in result['summaries'].items():
        lines.append(f'=== Logger {sn}: {len(S["files"])} log files ===')
        lines.append(f'Time: {S["start"]} to {S["end"]} UTC ({S["duration_hours"]:.2f} hours)')
        lines.append(f'Firmware: {", ".join(S["firmware"]) or "unknown"}; GPS fixes: {S["n_gps_fixes"]}')
        lines.append('Health summary (obviously bad values removed):')
        lines.append(S['stats'].to_string(float_format = lambda x: f'{x:.2f}', formatters = {'n': '{:.0f}'.format}))
        lines.append('')
        if len(S['warnings']) == 0:
            lines.append('No warnings.')
        else:
            lines.append(f'{len(S["warnings"])} warnings:')
            lines += [f'  WARNING: {w}' for w in S['warnings']]
        lines.append('')
    for p in result['problems']:
        lines.append(f'NOTE: {p}')
    return '\n'.join(lines)


def _moving_average(x, y, window = '30min'):
    """Centered moving average of y over a time window, ignoring NaN; x is a datetime Series."""
    s = pd.Series(y.values, index = pd.DatetimeIndex(x.values)).dropna().sort_index()
    m = s.rolling(window, center = True, min_periods = 1).mean()
    return m.index, m.values


def plot_aspen_summary(result, SN = None, filename = None):
    """
    Figure of temperature, relative humidity, battery voltage, and power draw (V * mA, in mW, log scale) against UTC
    time, with a 30-minute moving average on the power, from the output of summarize_aspen_logs. Obviously bad values are left
    out. Restarts show as vertical dotted lines.
    Health lines with no UTC time (before the first GPS fix of a logger that never got one) are not plotted.

    Parameters
    ----------
    result : dict
        Output of summarize_aspen_logs.
    SN : str, optional
        Serial number to plot (required if the result has more than one).
    filename : str, optional
        If given, the figure is saved here.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, NullFormatter
    if SN is None:
        if len(result['summaries']) != 1:
            raise ValueError('SN is required when there is more than one serial number')
        SN = next(iter(result['summaries']))
    H = result['summaries'][SN]['health']
    H = H[H['t'].notna()]
    if len(H) == 0:
        raise ValueError(f'{SN}: no health lines to plot')
    x = pd.to_datetime(H['t'], unit = 's')
    restarts = np.where(np.diff(H['tcxoSec'].values) < 0)[0] + 1
    H = H.assign(mW = H['V'] * H['mA'])
    fig, axes = plt.subplots(4, 1, sharex = True, figsize = (10, 8))
    panels = [(['temp'], 'Temperature (C)'), (['RH'], 'Relative humidity (%)'),
              (['batt', 'V'], 'Battery voltage (V)'), (['mW'], 'Power draw (mW)')]
    for ax, (cols, label) in zip(axes, panels):
        for col in cols:
            name = {'batt': 'batt (analog)', 'V': 'V (power monitor)'}.get(col, col)
            ax.plot(x, H[col], '.', ms = 1.5, alpha = 0.5, label = name)
            if col == 'mW':
                ax.plot(*_moving_average(x, H[col]), '-', lw = 1.5, color = 'k', label = '30-min mean')
        for i in restarts:
            if pd.notna(x.iloc[i]):
                ax.axvline(x.iloc[i], color = 'k', ls = ':', lw = 0.8)
        ax.set_ylabel(label)
        ax.legend(loc = 'best', markerscale = 6, ncol = 2, fontsize = 'small')
        if cols == ['mW']:
            ax.set_yscale('log')
            ## plain numbers (50, 100), not 5x10^1, at 1-2-3-5 steps
            lo, hi = ax.get_ylim()
            ticks = [m * 10.0 ** k for k in range(-3, 7) for m in (1, 2, 3, 5)]
            ax.set_yticks([t for t in ticks if lo <= t <= hi])
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'{v:g}'))
            ax.yaxis.set_minor_formatter(NullFormatter())
    axes[-1].set_xlabel('UTC time')
    axes[0].set_title(f'Logger {SN}: health summary')
    fig.autofmt_xdate()
    fig.tight_layout()
    if filename is not None:
        fig.savefig(filename, dpi = 120)
    return fig
