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
def check_lags2(DB, winlength = 1000, fl = 0.5, fh = 20, maxshift = 5):
    stations = DB.station.unique()
    from obspy.signal.cross_correlation import xcorr
    st = obspy.Stream()
    st.filter('bandpass', freqmin=fl, freqmax = fh)
    for fn in DB.filename:
        st += obspy.read(fn)
    t0 = min(DB.t1).replace(minute=0, second=0, microsecond=0)
    t = []
    lag = np.zeros([len(stations)-1, 1])
    xc_coef = np.zeros([len(stations)-1, 1])
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
            st_test = st.slice(t1, t1 + winlength)
            for i in range(len(stations) - 1):
                xc_output = xcorr(st_test.traces[0], st_test.traces[i+1], maxshift, full_xcorr=True)
                test_lags.append(xc_output[0])
                test_xc_coefs.append(xc_output[1])
            t.append(t1)
            lag = np.hstack([lag, np.array(test_lags).reshape(len(stations)-1, 1)])
            xc_coef = np.hstack([xc_coef, np.array(test_xc_coefs).reshape(len(stations)-1, 1)])
        except:
            pass
    return([t, lag, xc_coef])

def check_lags(DB, winlength = 1000, fl = 0.5, fh = 20, maxshift = 5):
    stations = DB.station.unique()
    from obspy.signal.cross_correlation import xcorr
    st = obspy.Stream()
    st.filter('bandpass', freqmin=fl, freqmax = fh)
    for fn in DB.filename:
        st += obspy.read(fn)
    t0 = min(DB.t1).replace(minute=0, second=0, microsecond=0)
    t = []
    lag = np.zeros([len(stations)-1, 1])
    xc_coef = np.zeros([len(stations)-1, 1])
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
            st.trim(t1)
            st_test = st.copy()
            st_test.trim(t1, t1 + winlength)
            for i in range(len(stations) - 1):
                xc_output = xcorr(st_test.traces[0], st_test.traces[i+1], maxshift, full_xcorr=True)
                test_lags.append(xc_output[0])
                test_xc_coefs.append(xc_output[1])
            t.append(t1)
            lag = np.hstack([lag, np.array(test_lags).reshape(len(stations)-1, 1)])
            xc_coef = np.hstack([xc_coef, np.array(test_xc_coefs).reshape(len(stations)-1, 1)])
        except:
            pass
    return([t, lag, xc_coef])

    
## run from ~/Work/Gem_Tests/2020-09-30_BalconyTests
#import CheckDiscontinuity
import gemlog
#DB = gemlog.make_db('mseed', '2020-09-30T0[3-9]*{058,061,077}*')
DB = gemlog.make_db('mseed', '*')
import glob
#CheckDiscontinuity(glob.glob('mseed/'))
lag_output = check_lags(DB)
#lag_output = gemlog_aux.lag_output
import matplotlib.pyplot as plt
N = lag_output[1].shape[0]
plt.subplot(2,1,1)
for i in range(N):
    plt.plot(lag_output[1][i,:]+i/6, '.')

plt.vlines(np.arange(14) * 84.6, -5, 5)
plt.subplot(2,1,2)
for i in range(N):
    plt.plot(lag_output[2][i,:], '.')

plt.vlines(np.arange(14) * 84.6, 0, 1)
plt.show()
import pickle
with open('lag_output_1000.pkl', 'wb') as file:
    pickle.dump(lag_output, file)
    
    
#fn = glob.glob('mseed/*095*')
#CheckDiscontinuity(fn)



#import sys
#sys.path.append('/home/jake/Dropbox/Gem_logger/gem_package/python/gemlog_python/gemlog/gemlog')
#import gemlog
#gemlog.convert(rawpath = 'raw', SN = '095')
#fn = glob.glob('raw/*')
#fn.sort()
#L=gemlog.ReadGem_v0_9(fn)
#L2=gemlog.ReadGem(path='raw', SN ='095')


