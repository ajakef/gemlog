import numpy as np
import matplotlib.pyplot as plt
import obspy
from gemlog import *
def PlotAmp(DB):
    allSta = DB.station.unique()
    allSta.sort()
    for sta in DB.station.unique():
        w = np.where(DB.station == sta)[0]
        w.sort()
        plt.plot(DB.t1.iloc[w], np.log10(DB.amp_HP.iloc[w]), '.')
        print(str(sta) + ' ' + str(np.quantile(DB.amp_HP.iloc[w], 0.25)))
    plt.legend(allSta)
    plt.show()

import obspy
def CheckDiscontinuity(files):
    st = obspy.Stream()
    for fn in files:
        st += obspy.read(fn)
    st.merge()
    st = st.split()
    if len(st) > 1:
        for tr in st:
            print(tr.id + ' ' + tr.stats.starttime.isoformat() + '--' + tr.stats.endtime.isoformat())

import gemlog, obspy
import numpy as np

def check_lags(DB, winlength = 1000, fl = 0.5, fh = 20, maxshift = 10):
    stations = DB.station.unique()
    nsta = len(stations)
    from obspy.signal.cross_correlation import xcorr, correlate, xcorr_max
    st = obspy.Stream()
    st.filter('bandpass', freqmin=fl, freqmax = fh)
    for fn in DB.filename:
        st += obspy.read(fn)
    t0 = min(DB.t1).replace(minute=0, second=0, microsecond=0)
    t = []
    lag = np.zeros([nsta, 1])
    xc_coef = np.zeros([nsta, 1])
    consistency = []
    N = int((max(DB.t2) - t0)/winlength)
    count = 0
    t1 = t0 + 1e-6
    #while t1 < (t0 + 10*winlength):#max(DB.t2):
    while t1 < max(DB.t2):
        count += 1
        t1 += winlength
        print(str(count) + ' of ' + str(N))
        try:
            test_lags = []
            test_xc_coefs = []
            #st.trim(t1)
            st_test = st.slice(t1, t1 + winlength)
            for i in range(nsta):
                xc_output = xcorr_max(correlate(st_test.traces[i], st_test.traces[(i+1) % nsta], maxshift))
                test_lags.append(xc_output[0])
                test_xc_coefs.append(xc_output[1])
            t.append(t1)
            lag = np.hstack([lag, np.array(test_lags).reshape(nsta, 1)])
            xc_coef = np.hstack([xc_coef, np.array(test_xc_coefs).reshape(nsta, 1)])
            consistency.append(np.sum(test_lags))
        except:
            pass
    return([t, lag[:,1:], xc_coef[:,1:], np.array(consistency)])


def plot_lags(lag0, lag1, use_consistency = True):
    import matplotlib.pyplot as plt
    N=3
    ## Plot lags for the two conversions
    plt.subplot(2,2,1)
    ## consistency multiplier: 1 if consistent, NaN otherwise
    if use_consistency:
        c1 = (np.array(lag1[3])==0) / (np.array(lag1[3])==0)
        c0 = (np.array(lag0[3])==0) / (np.array(lag0[3])==0)
    else:
        c1 = 1
        c0 = 1
    for i in range(N):
        plt.plot(c1 * lag1[1][i,:]+i/6, '.')
    
    plt.subplot(2,2,3)
    for i in range(N):
        plt.plot(c0 * lag0[1][i,:]+i/6, '.')

    ## Plot the correlation coefficient
    plt.subplot(2,2,2)
    for i in range(N):
        plt.plot(c1*lag1[2][i,:], '.')

    plt.subplot(2,2,4)
    for i in range(N):
        plt.plot(c0 * lag0[2][i,:], '.')

    plt.legend([str(i) for i in range(N)], loc = 'lower right')
    plt.show()
