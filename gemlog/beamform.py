"""
Fast conventional beamforming, a drop-in replacement for
obspy.signal.array_analysis.array_processing.

obspy's version takes time proportional to (number of stations)^2 for
beamforming (method=0) and builds a very large table of steering vectors. This
version takes time proportional to the number of stations, with the same
inputs and outputs. Capon (method=1) is passed to obspy unchanged.
"""
import math

import numpy as np
from matplotlib.dates import datestr2num
from obspy.signal.array_analysis import (
    array_processing as _obspy_array_processing,
    get_geometry,
    get_spoint,
)
from obspy.signal.invsim import cosine_taper
from obspy.signal.util import next_pow_2


def _beamformer_bf_setup(geometry, sll_x, sll_y, sl_s, grdpts_x, grdpts_y, freqs):
    # The time shift for station i is sx*x_i + sy*y_i, so the phase factor is a
    # product of one factor depending on sx and one depending on sy.
    sx = sll_x + np.arange(grdpts_x) * sl_s
    sy = sll_y + np.arange(grdpts_y) * sl_s
    # single precision, to match obspy's time shift table
    tx = np.require(np.outer(geometry[:, 0], sx), dtype=np.float32)
    ty = np.require(np.outer(geometry[:, 1], sy), dtype=np.float32)
    phase = 2j * np.pi * np.asarray(freqs, dtype=np.float64)[:, None, None]
    return np.exp(phase * tx[None]), np.exp(phase * ty[None])


def _beamformer_bf(ft, phase_x, phase_y, nstat, prewhiten, dpow, relpow_map, abspow_map):
    # P(f) = e.H R e with R = ft ft.H has rank one, so P(f) = |e.H ft|^2
    relpow_map.fill(0.)
    abspow_map.fill(0.)
    nf = ft.shape[1]
    for n in range(nf):
        beam = (phase_x[n] * ft[:, n, np.newaxis]).T @ phase_y[n]
        pow_n = beam.real ** 2 + beam.imag ** 2
        abspow_map += pow_n
        if prewhiten == 1:
            relpow_map += pow_n / (pow_n.max() * nf * nstat)
    if prewhiten != 1:
        relpow_map += abspow_map / dpow


def array_processing(stream, win_len, win_frac, sll_x, slm_x, sll_y, slm_y,
                     sl_s, semb_thres, vel_thres, frqlow, frqhigh, stime,
                     etime, prewhiten, verbose=False, coordsys='lonlat',
                     timestamp='mlabday', method=0, store=None):
    """
    Same arguments and output as obspy.signal.array_analysis.array_processing;
    see the obspy documentation. Much faster for method=0 (beamforming) with
    many stations.
    """
    if method != 0:
        return _obspy_array_processing(
            stream, win_len, win_frac, sll_x, slm_x, sll_y, slm_y, sl_s,
            semb_thres, vel_thres, frqlow, frqhigh, stime, etime, prewhiten,
            verbose=verbose, coordsys=coordsys, timestamp=timestamp,
            method=method, store=store)
    res = []
    eotr = True
    fs = stream[0].stats.sampling_rate
    if len(stream) != len(stream.select(sampling_rate=fs)):
        raise ValueError('in sonic sampling rates of traces in stream are not equal')
    grdpts_x = int(((slm_x - sll_x) / sl_s + 0.5) + 1)
    grdpts_y = int(((slm_y - sll_y) / sl_s + 0.5) + 1)
    geometry = get_geometry(stream, coordsys=coordsys, verbose=verbose)
    spoint, _epoint = get_spoint(stream, stime, etime)
    nstat = len(stream)
    nsamp = int(win_len * fs)
    nstep = int(nsamp * win_frac)
    nfft = next_pow_2(nsamp)
    deltaf = fs / float(nfft)
    nlow = max(1, int(frqlow / float(deltaf) + 0.5))  # avoid the offset
    nhigh = min(nfft // 2 - 1, int(frqhigh / float(deltaf) + 0.5))  # avoid nyquist
    nf = nhigh - nlow + 1
    phase_x, phase_y = _beamformer_bf_setup(
        geometry, sll_x, sll_y, sl_s, grdpts_x, grdpts_y,
        (nlow + np.arange(nf)) * deltaf)
    ft = np.empty((nstat, nf), dtype=np.complex128)
    newstart = stime
    tap = cosine_taper(nsamp, p=0.22)
    offset = 0
    relpow_map = np.empty((grdpts_x, grdpts_y), dtype=np.float64)
    abspow_map = np.empty((grdpts_x, grdpts_y), dtype=np.float64)
    while eotr:
        try:
            for i, tr in enumerate(stream):
                dat = tr.data[spoint[i] + offset: spoint[i] + offset + nsamp]
                dat = (dat - dat.mean()) * tap
                ft[i, :] = np.fft.rfft(dat, nfft)[nlow:nlow + nf]
        except IndexError:
            break
        dpow = nstat * np.abs(ft * ft.conj()).sum()
        _beamformer_bf(ft, phase_x, phase_y, nstat, prewhiten, dpow,
                       relpow_map, abspow_map)
        ix, iy = np.unravel_index(relpow_map.argmax(), relpow_map.shape)
        relpow, abspow = relpow_map[ix, iy], abspow_map[ix, iy]
        if store is not None:
            store(relpow_map, abspow_map, offset)
        slow_x = sll_x + ix * sl_s
        slow_y = sll_y + iy * sl_s
        slow = max(np.sqrt(slow_x ** 2 + slow_y ** 2), 1e-8)
        azimut = 180 * math.atan2(slow_x, slow_y) / math.pi
        baz = azimut % -360 + 180
        if relpow > semb_thres and 1. / slow > vel_thres:
            res.append(np.array([newstart.timestamp, relpow, abspow, baz, slow]))
            if verbose:
                print(newstart, (newstart + (nsamp / fs)), res[-1][1:])
        if (newstart + (nsamp + nstep) / fs) > etime:
            eotr = False
        offset += nstep
        newstart += nstep / fs
    res = np.array(res)
    if timestamp == 'julsec':
        pass
    elif timestamp == 'mlabday':
        res[:, 0] = res[:, 0] / (24. * 3600) + datestr2num('1970-01-01')
    else:
        raise ValueError("Option timestamp must be one of 'julsec', or 'mlabday'")
    return res
