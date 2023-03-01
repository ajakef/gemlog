from obspy.signal.cross_correlation import correlate, xcorr_max
import numpy as np
import pandas as pd
import glob, os
from gemlog.gem_network import _unique

xcorr_all(pattern = '*00..1*[0-3]', t1 = '2020-04-20', t2 = '2020-04-22')

st = obspy.read('2020-04-20*00..10[0-3]*')
output = xcorr_one_day(st)

def xcorr_all(path = '.', pattern = '*', t1 = '1970-01-01', t2 = '9999-12-31', IDs = None):
    return loop_through_days(xcorr_one_day, path, pattern, t1, t2, IDs)
    


def loop_through_days(function, path = '.', pattern = '*', t1 = '1970-01-01', t2 = '9999-12-31', IDs = None):
    try:
        t1 = obspy.UTCDateTime(t1)
    except:
        raise Exception('invalid t1')
    try:
        t2 = obspy.UTCDateTime(t2)
    except:
        raise Exception('invalid t2')
    
    ## make a database of files
    file_metadata = {'filename':[], 't1':[], 't2':[], 'ID':[]}
    filenames = glob.glob(os.path.join(path, pattern))
    for filename in filenames:
        st = obspy.read(filename, headonly=True)
        for tr in st:
            file_metadata['filename'].append(filename)
            file_metadata['t1'].append(tr.stats.starttime)
            file_metadata['t2'].append(tr.stats.endtime)
            file_metadata['ID'].append(tr.id)

    file_metadata_df = pd.DataFrame.from_dict(file_metadata)

    t1 = np.max([t1, file_metadata_df.t1.min()])
    t2 = np.min([t2, file_metadata_df.t2.max()])
    if IDs is None:
        IDs = _unique(file_metadata_df.ID)

    rows_to_keep = np.where((file_metadata_df.t2 >= t1) & \
                    (file_metadata_df.t1 <= t2) & \
                    np.isin(file_metadata_df.ID, IDs))[0]
    file_metadata_df = file_metadata_df.iloc[rows_to_keep,:]
    file_metadata_df.sort_values('t1', ignore_index=True, inplace=True)

    ## loop through days. be careful to avoid funny business with leap seconds.
    day_starts = [t1]
    while True:
        new_time = (day_starts[-1]+86400+2).replace(hour=0, minute=0, second=0) # add 24 hours + 2 sec, then round down
        if new_time < t2:
            day_starts.append(new_time)
        else:
            break
    day_starts = np.array(day_starts)

    day_ends = [day_starts[1]]
    while True:
        new_time = (day_ends[-1]+86400+2).replace(hour=0, minute=0, second=0) # add 24 hours + 2 sec, then round down
        if new_time < t2:
            day_ends.append(new_time)
        else:
            day_ends.append(t2)
            break
    day_ends = np.array(day_ends)

    files_read = []
    output_list = []
    st = obspy.Stream()
    for (day_start, day_end) in zip(day_starts, day_ends):
        indices = np.where((file_metadata_df.t1 <= day_end) &
                           (file_metadata_df.t2 >= day_start))[0]
        ## for efficiency, avoid reading files twice if they cover multiple days or IDs
        for filename in file_metadata_df.filename[indices]:
            if filename not in files_read:
                st += obspy.read(filename)
                files_read.append(filename)

        st.merge()
        ## throw out data before the start of this day, and any unused traces
        st.trim(day_start, t2)
        for tr in st:
            if tr.id not in IDs:
                st.pop(tr)

        ## finally, apply whatever function you have to the data
        ## the function will return a pd.DataFrame. append to a list, then merge at the end
        day_output = function(st.slice(day_start, day_end))
        output_list.append(day_output)
    ## done looping. merge the output and return.
    return pd.concat(output_list, ignore_index = True)
        
        
def xcorr_one_day(st, fl = 1, fh = 40, win_length_sec = 10, overlap = 0):
    st.detrend('linear')
    st.filter('bandpass', freqmin = fl, freqmax = fh)
    output = apply_function_windows(st, xcorr_function, win_length_sec, overlap)

    ## reformat UTCDateTimes as string
    #output['t_mid'] = [t.isoformat() for t in output['t_mid']]
    return pd.DataFrame.from_dict(output)

def xcorr_function(st, maxshift_seconds = 1):
    print(st)
    dt = st[0].stats.delta
    st.detrend('linear')
    st.taper(0.05, 'hann') # Hann window, default taper for SAC
    output_dict = {'mean_coef':0}
    consistency = 0
    ## simplest loop: go through sensors in order, O(N). Not robust if a single sensor misbehaves.
    for i in range(0, len(st)):
        j = (i+1) % len(st)
        tr1 = st[i]
        tr2 = st[j]
        output_dict[f'rms_{tr1.stats.station}_{tr1.stats.location}'] = tr1.std()
        pair_name = f'{tr1.stats.station + tr1.stats.location}_{tr2.stats.station + tr2.stats.location}'
        shift, value = xcorr_max(correlate(tr1.data, tr2.data, int(np.round(maxshift_seconds / dt))), abs_max = False)
        output_dict[f'lag_{pair_name}'] = shift * dt
        output_dict[f'r_{pair_name}'] = value
        consistency += shift/len(st)
        output_dict['mean_coef'] += value/len(st)
    ## consistency: allow up to one sample error per cross-correlation
    output_dict['consistency'] = np.abs(consistency) <= len(st)
    return output_dict



def apply_function_windows(st, f, win_length_sec, overlap = 0.5):
    """
    Run an analysis (or suite of analyses) on overlapping windows for some dataset
    
    Parameters:
    -----------
    st : obspy.Stream
    Stream including data to divide into windows and analyze

    f : function
    Accepts single variable "st" (obspy.Stream), returns dictionary of results

    win_length_sec : float
    Length of analysis windows [seconds]

    overlap : float
    Proportion of a window that overlaps with the previous window [unitless, 0-1]
 
    Returns:
    --------
    dictionary with following items:
    --t_mid (obspy.UTCDateTime): mean time of each window
    --all elements of the output of "f", joined into numpy arrays

    Note:
    -----
    For each time window, the stream is trimmed to fit, but not tapered, detrended, or otherwise 
    processed. If those steps are necessary, be sure they are carried out in f().
"""
    # f must input a stream and return a dict
    eps = 1e-6
    t1 = st[0].stats.starttime
    t2 = st[0].stats.endtime
    data_length_sec = t2 - t1
    num_windows = 1 + int(np.ceil((data_length_sec - win_length_sec) / (win_length_sec * (1 - overlap)) - eps))
    output_dict = {'t_mid':[]}
    print(num_windows)
    for i in range(num_windows):
        win_start = t1 + i*(data_length_sec - win_length_sec) / (num_windows-1)
        st_tmp = st.slice(win_start-eps, win_start + win_length_sec - eps, nearest_sample = False)
        win_dict = f(st_tmp)
        #if i == 0:
        #    output_dict = {key:[] for key in win_dict.keys()}
        #    output_dict['t_mid'] = []
        #for key in win_dict.keys():
        #    output_dict[key].append(win_dict[key])

        ## append data from this window to the output
        for key in win_dict.keys():
            if key not in output_dict.keys():
                output_dict[key] = list(np.repeat(np.nan, i))
            output_dict[key].append(win_dict[key])
        ## if a field from the output is missing in this window's data, append nan as a placeholder
        for key in output_dict.keys():
            if key not in win_dict.keys() and key != 't_mid':
                output_dict[key].append(np.nan)

            
        output_dict['t_mid'].append(win_start + win_length_sec/2)
    output_dict = {key:np.array(output_dict[key]) for key in output_dict.keys()}
    return output_dict
