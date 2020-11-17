import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import obspy, glob, gemlog
#from gemlog import *
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


def check_lags(DB, winlength = 1000, fl = 0.5, fh = 20, maxshift = 10, verbose = False):
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
        if verbose: print(str(count) + ' of ' + str(N))
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

def make_db(path, pattern = '*', savefile = None, verbose = False):
    """Create a database summarizing a set of converted data files.

    Parameters
    ----------
    path : str
        Path to folder containing converted data files to summarize.
    
    pattern : str, default '*'
        Glob-type pattern for selecting converted data files to summarize

    savefile : str, default None
        File name where database is written. Use 'savefile = None' to not save an output file.

    verbose : bool, default False
        If True, print progress updates.

    Returns
    -------
    pandas.DataFrame containing converted file database.
    """
    #path = 'mseed'
    #pattern = '*'
    files = glob.glob(path + '/' + pattern)
    files.sort()
    DB = []
    count = 0
    for file in files:
        tr = obspy.read(file)[0]
        maxVal = tr.data.max()
        minVal = tr.data.min()
        tr.detrend('linear')
        tr.filter('highpass', freq=0.5)
        amp_HP = tr.std()
        row = pd.DataFrame([[file, tr.stats.station, tr.stats.location, amp_HP, maxVal, minVal, tr.stats.starttime, tr.stats.endtime]], columns = ['filename', 'station', 'location', 'amp_HP', 'max', 'min', 't1', 't2'])
        DB.append(row)
        if((count % 100) == 0 and verbose):
            print(str(count) + ' of ' + str(len(files)))
        count = count + 1
    DB = pd.concat(DB)
    if savefile is not None:
        DB.to_csv(savefile)
    return(DB)

################################################

################

def calc_channel_stats(DB, t1, t2):
    """
    Calculate uptime and other statistics for all channels in a database.

    Parameters
    ----------
    DB : pandas.DataFrame
        Output of make_db().
    
    t1 : time-like 
        Start time for which statistics should be calculated.

    t2 : time-like 
        End time for which statistics should be calculated.

    Returns
    -------
    pandas.DataFrame containing the following columns:

        - station : station name
        - goodData : proportion of time window (t1-t2) that is not obviously bad (e.g., clipped)
        - anyData : proportion of time window (t1-t2) for which data are available
        - q1 : first quartile amplitude 
        - q3 : third quartile amplitude
    """
    import obspy, glob
    import pandas as pd
    from obspy import UTCDateTime as T
    #t1 = '2020-04-14'
    #t2 = '2020-04-24T20:00:00'
    t1 = obspy.core.UTCDateTime(t1)
    t2 = obspy.core.UTCDateTime(t2)
    numHour = (t2 - t1)/3600.0
    DB.t1 = DB.t1.apply(T)
    DB.t2 = DB.t2.apply(T)
    DB.goodData = (DB.amp_HP > 0.5) & (DB.amp_HP < 2e4) & ((DB.t2 - DB.t1) > 3598) & ((DB.t2 - DB.t1) < 3602)
    DB.anyData = (DB.amp_HP > 0) 
    out = []
    for sta in DB.station.unique():
        w = np.where((DB.station == sta) & (DB.t1 > t1) & (DB.t2 < t2))[0]
        if(len(w) == 0):
            continue
        else:
            q1 = np.quantile(np.array(DB.amp_HP)[w], 0.25)
            q3 = np.quantile(np.array(DB.amp_HP)[w], 0.75)
            out.append(pd.DataFrame([[sta, np.sum(np.array(DB.goodData)[w])/numHour, np.sum(np.array(DB.anyData)[w])/numHour, q1, q3]], columns = ['station', 'goodData', 'anyData', 'q1', 'q3']))
    out = pd.concat(out)
    return(out)

