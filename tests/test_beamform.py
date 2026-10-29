import numpy as np
import obspy
import pytest
from obspy.signal.array_analysis import array_processing as obspy_array_processing

from gemlog.beamform import array_processing


def make_stream(nstat, coordsys, seed):
    rng = np.random.RandomState(seed)
    st = obspy.Stream()
    for _ in range(nstat):
        tr = obspy.Trace(rng.randn(50 * 40))
        tr.stats.sampling_rate = 50
        if coordsys == 'xy':
            coords = dict(x=rng.uniform(-0.5, 0.5), y=rng.uniform(-0.5, 0.5), elevation=0.)
        else:
            coords = dict(longitude=rng.uniform(0, 0.01), latitude=rng.uniform(0, 0.01), elevation=0.)
        tr.stats.coordinates = obspy.core.AttribDict(coords)
        st += tr
    return st


@pytest.mark.parametrize('coordsys', ['xy', 'lonlat'])
@pytest.mark.parametrize('prewhiten', [0, 1])
@pytest.mark.parametrize('nstat', [3, 12])
def test_matches_obspy(nstat, prewhiten, coordsys):
    st = make_stream(nstat, coordsys, nstat)
    args = (st, 4, 0.5, -2, 2, -1, 3, 0.25, -1e99, -1e99, 1, 10,
            st[0].stats.starttime, st[0].stats.endtime, prewhiten)
    for timestamp in ('mlabday', 'julsec'):
        ref = obspy_array_processing(*args, coordsys=coordsys, timestamp=timestamp)
        out = array_processing(*args, coordsys=coordsys, timestamp=timestamp)
        assert out.shape == ref.shape
        np.testing.assert_allclose(out, ref, rtol=1e-5)


def test_thresholds_and_capon_passthrough():
    st = make_stream(5, 'xy', 0)
    args = (st, 4, 0.5, -2, 2, -2, 2, 0.25, 0.01, 0.1, 1, 10,
            st[0].stats.starttime, st[0].stats.endtime, 0)
    np.testing.assert_allclose(array_processing(*args, coordsys='xy'),
                               obspy_array_processing(*args, coordsys='xy'), rtol=1e-5)
    np.testing.assert_allclose(array_processing(*args, coordsys='xy', method=1),
                               obspy_array_processing(*args, coordsys='xy', method=1))
