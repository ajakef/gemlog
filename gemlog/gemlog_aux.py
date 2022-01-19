import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import obspy, glob, gemlog, os
import scipy.integrate
import scipy.interpolate
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
        try:
            tr = obspy.read(file)[0]
        except:
            print('failed to read ' + file + ', skipping')
            continue
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

_noise_path = os.path.join(os.path.dirname(__file__), "data", "noise")

def gem_noise(freq = None, spectype = 'power', version = '1.0', freq_min = None, freq_max = 50):
    """Calculate Gem self-noise as a one-sided spectrum, and integrated over 
    a defined frequency band

    Parameters
    ----------
    freq : numpy array
        Input frequencies for spectrum. If None, use default.
    spectype : str
        Type of spectrum and noise calculations: either 'power', 'amp', or 'dB' 
    version : str
        Gem version to find noise spectrum for
    freq_min : float
        minimum frequency of band for calculating noise
    freq_max : float
        maximum frequency of band for calculating noise
    
    Returns
    -------
    Dictionary with following items:
    freqs: frequencies of output spectrum
    spectrum: one-sided self-noise spectrum
    type: type of spectrum (power, amplitude, or dB, depends on spectrype input)
    spectrum_units: units of spectrum (depends on spectype input)
    freq_min: lower bound of noise band
    freq_max: upper bound of noise band
    noise: integrated noise over band
    noise_units: untis of integrated noise (depends on spectype input)

    Example
    -------
    ## calculate gem noise over the 0.1-20 Hz band as an RMS amplitude
    noise_info = gem_noise(version = '1.0', freq_min = 0.1, freq_max = 20, spectype = 'amp')
    print(noise_info['noise'])
    print(noise_info['noise_units'])
    
    """
    # use cases:
    # 1: wants output spectrum for given frequency vector
    # 2: wants output frequencies and spectrum with no frequency input
    # in order to permit 2, we need to return frequencies and spectrum
    if(float(version) >= 0.97):
        noise_spec = pd.read_csv(os.path.join(_noise_path, 'Gem_v0.98_Noise_spec.txt'))
    else:
        raise Exception('version %s not supported' % version)
    freq_in = noise_spec['f']
    power_in = noise_spec['amp']**2

    return _noise_spectrum_helper(freq_in, power_in, freq, spectype, freq_min, freq_max)

def ims_noise(model = 'low', freq = None, spectype = 'power', freq_min = None, freq_max = 50):
    """Return noise from an IMS noise model as a one-sided spectrum integrated
    over a defined frequency band

    Parameters
    ----------
    freq : numpy array
        Input frequencies for spectrum. If None, use default.
    spectype : str
        Type of spectrum and noise calculations: either 'power', 'amp', or 'dB' 
    freq_min : float
        minimum frequency of band for calculating noise
    freq_max : float
        maximum frequency of band for calculating noise
    
    Returns
    -------
    Dictionary with following items:
    freqs: frequencies of output spectrum
    spectrum: low noise model (one-sided self-noise spectrum)
    type: type of spectrum (power, amplitude, or dB, depends on spectrype input)
    spectrum_units: units of spectrum (depends on spectype input)
    freq_min: lower bound of noise band
    freq_max: upper bound of noise band
    noise: integrated noise over band
    noise_units: units of integrated noise (depends on spectype input)

    Example
    -------
    ## calculate gem noise over the 0.1-20 Hz band as an RMS amplitude
    noise_info = ims_noise(model = 'low', freq_min = 0.1, freq_max = 20, spectype = 'amp')
    print(noise_info['noise'])
    print(noise_info['noise_units'])
    
    """
    # use cases:
    # 1: wants output spectrum for given frequency vector
    # 2: wants output frequencies and spectrum with no frequency input
    # in order to permit 2, we need to return frequencies and spectrum
    noise_spec = pd.read_csv(os.path.join(_noise_path, 'IMSNOISE_MIN_MED_MAX.txt'),
                             names = ['f', 'min', 'med', 'max'], sep = r'\s+')
    freq_in = noise_spec['f']
    if model.lower() in ['low', 'min']:
        power_in = 10**noise_spec['min']
    elif model.lower() == 'med':
        power_in = 10**noise_spec['med']
    elif model.lower() in ['high', 'hi', 'max']:
        power_in = 10**noise_spec['max']
    else:
        raise Exception('Invalid IMS model: %s' % model)
    
    return _noise_spectrum_helper(freq_in, power_in, freq, spectype, freq_min, freq_max)



def _noise_spectrum_helper(freq_in, power_in, freq_out, spectype, freq_min, freq_max):
    if freq_out is None:
        freq_out = freq_in
    spec_function = scipy.interpolate.CubicSpline(freq_in, power_in, extrapolate = False)
    power = spec_function(freq_out)
    if (freq_min is not None) and (freq_max is not None):
        noise = scipy.integrate.quad(spec_function, freq_min, freq_max)[0] # 'quad' outputs are the estimated integral and error bar
    else:
        noise = np.NaN

    ## decide what kind of spectrum to return (default "power")
    if spectype.lower() == 'power':
        spec = power
        spec_units = 'Pa^2/Hz'
        noise_units = 'Pa^2'
    elif spectype.lower()[:3] == 'amp':
        spec = np.sqrt(power)
        noise = np.sqrt(noise)
        spec_units = 'Pa/sqrt(Hz)'
        noise_units = 'Pa'
    elif spectype.lower() == 'db':
        spec = 10*np.log10(power)
        noise = 10*np.log10(noise)
        spec_units = 'dB Pa/sqrt(Hz)'
        noise_units = 'dB Pa'
    else:
        raise Exception('spectype %s not supported' % spectype)

    return {'freqs':freq_out, 'spectrum':spec, 'type':spectype, 'spectrum_units':spec_units, 'freq_min':freq_min, 'freq_max':freq_max, 'noise':noise, 'noise_units':noise_units}


def _interpolate_stream(st, gap_limit_sec = 0.1):
    ## look for short data gaps and interpolate through them
    ## st must consist of data from just a single station
    st.merge()
    st = st.split()

    ## If there are no data gaps, the stream will only contain a single trace. Don't change anything.
    if len(st) == 1:
        return st
    
    for i in range(len(st) - 1):
        t1 = st[i].stats.endtime
        t2 = st[i+1].stats.starttime
        if (t2 - t1) < gap_limit_sec:
            ## define a new trace consisting of a masked array of the data gap and immediate surroundings
            st_new = st.slice(t1 - gap_limit_sec, t2 + gap_limit_sec)
            st_new.merge()
            tr_new = st_new[0]

            ## interpolate the missing values using a cubic spline
            data_good = ~tr_new.data.mask
            t = np.arange(len(tr_new.data))
            interp_function = scipy.interpolate.CubicSpline(t[data_good], tr_new.data[data_good], extrapolate = False)
            tr_new.data[~data_good] = interp_function(t[~data_good])

            ## add the new trace back to the stream
            st = st + tr_new
        else:
            ## if the data gap is not short, it will not be interpolated. For now, do nothing.
            pass
        
    st.merge() ## merge all the traces (done in place)
    st = st.split() ## if there are still data gaps, break them up (ensuring that the result is not a masked array)

    return st

            
            
            
    
