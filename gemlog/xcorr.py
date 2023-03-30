from obspy.signal.cross_correlation import correlate, xcorr_max
import obspy
import numpy as np
import pandas as pd
import glob, os, traceback, sys, getopt, argparse, re
from gemlog.gem_network import _unique
import gemlog
import scipy.interpolate

def xcorr_all_terminal(input = sys.argv[1:]):
    examples_text = f'''
    Examples: (replace '/' with '\' if using Windows)\n
    # process all data in mseed_data between 2022-09-01 00:00:00 UTC and 05:00:00 UTC from stations 121, 122, and 123
    waveform_calc_lags -1 2022-09-01 -2 2022-09-01_05:00:00 -i 121,122,123 -o output_file.csv mseed_data/* 

    # process all data in mseed_data except stations 100 and 110, first filtering between 10-20 Hz instead of default frequencies
    waveform_calc_lags -x 100,110 -L 10 -H 20 -o 10-20Hz_output_file.csv mseed_data/* 

    # process all data in mseed_data, upsampling by 4x to improve precision, and using a 30-second window length
    waveform_calc_lags -u 4 -w 30 -o upsampled_30sec_output_file.csv mseed_data/*    

    gemlog version {gemlog.__version__}'''
    parser = argparse.ArgumentParser(description='Use cross-correlation to find delays between waveform files, and calculate backazimuth and horizontal slowness.',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog = examples_text)
    parser.add_argument('files', nargs='+', help='List of data files to process (wildcards are allowed)')
    
    parser.add_argument('-o', '--output_file', nargs = 1, default=None, help='Output file to write')
    parser.add_argument('-i', '--include_IDs', help='Station IDs to include in processing (default: all)')
    parser.add_argument('-x', '--exclude_IDs', help='Station IDs to exclude from processing (default: none)')
    parser.add_argument('-1', '--t1', default='1970-01-01', help='Time to start processing data (default: beginning of data)')
    parser.add_argument('-2', '--t2', default='9999-12-31', help='Time to stop processing data (default: end of data)')
    parser.add_argument('-L', '--freq_low', default=5, help='Low corner frequency (default: 5 Hz)', type = float)
    parser.add_argument('-H', '--freq_high', default=40, help='High corner frequency (default: 40 Hz)', type = float)
    parser.add_argument('-w', '--window_length_seconds', default=5, help='Length of time windows to cross-correlate (default: 5 seconds)', type = float)
    parser.add_argument('-p', '--overlap', default=0, help='Window overlap fraction; 0 <= overlap < 1 (default: 0)', type = float)
    parser.add_argument('-u', '--upsample_ratio', default=1, help='Ratio by which waveform data should be upsampled to improve time resolution (default: 1 for no upsampling)', type = float)
    parser.add_argument('-q', '--quiet', action = 'store_true')
    args = parser.parse_args(input)

    if args.include_IDs is not None and type(args.include_IDs) is str:
        include_IDs = args.include_IDs.split(',')
    else:
        include_IDs = None
    if (args.exclude_IDs is not None) and type(args.exclude_IDs) is str:
        exclude_IDs = args.exclude_IDs.split(',')
    else:
        exclude_IDs = None

    if args.output_file is None:
        raise(Exception('output_file is a required input'))

    try:
        obspy.UTCDateTime(args.t1)
    except:
        raise Exception(f'Start time {args.t1} could not be interpreted as a UTC date-time')
    try:
        obspy.UTCDateTime(args.t2)
    except:
        raise Exception(f'End time {args.t2} could not be interpreted as a UTC date-time')


    ## Print the arguments. For conciseness, skip input files. To avoid confusion, skip any "None".
    print('')
    args_dict = vars(args).copy()
    args_dict.pop('files')
    args_to_print = {key:value for key, value in args_dict.items() if value is not None}
    print(args_to_print)
    print('')
    
    xcorr_df = xcorr_all(args.files, t1 = args.t1, t2 = args.t2, IDs = include_IDs,
                         exclude_IDs = exclude_IDs, fl = args.freq_low, fh = args.freq_high,
                         win_length_sec = args.window_length_seconds, overlap = args.overlap,
                         upsample_ratio = args.upsample_ratio, quiet = args.quiet)
    xcorr_df.to_csv(args.output_file[0], sep = ',', index = False)

    

def calculate_direction_terminal(input = sys.argv[1:]):
    epilog = f'gemlog version {gemlog.__version__}'
    parser = argparse.ArgumentParser(description='Invert time lags found by cross-correlation to find backazimuth and horizontal slowness.', formatter_class=argparse.HelpFormatter, epilog = epilog)
    parser.add_argument('-i', '--input_file', nargs = 1, help='File with time lags created by waveform_xc')
    parser.add_argument('-l', '--location_file', nargs = 1, help='stationXML file containing station locations')
    parser.add_argument('-o', '--output_file', nargs = 1, help='Output file to write, including azimuth and horizontal slowness results')
    args = parser.parse_args(input)

    print(args)
    

    if args.input_file is None:
        raise(Exception('input_file is a required input'))
    if args.location_file is None:
        raise(Exception('location_file is a required input'))
    if args.output_file is None:
        raise(Exception('output_file is a required input'))
    
    ## check that the input files exist and are readable
    try:
        xcorr_df = pd.read_csv(args.input_file[0])
    except:
        raise Exception(f'Failed to open input_file {args.input_file[0]}')
    try:
        locations = get_coordinates(obspy.read_inventory(args.location_file[0]))
    except:
        raise Exception(f'Failed to open location_file {args.location_file[0]}')

    ## check that the output file is writeable (do this before the inversion in case of long runtime)
    try:
        with open(args.output_file[0], 'w') as f:
            pass
    except:
        raise Exception(f'Cannot write to output_file {args.output_file[0]}')

    results = invert_for_slowness(xcorr_df, locations)
    results.to_csv(args.output_file[0], sep = ',', index = False)

    
def xcorr_all(files, t1 = '1970-01-01', t2 = '9999-12-31', IDs = None, exclude_IDs = None, fl = 5, fh = 20, win_length_sec = 5, overlap = 0, upsample_ratio = 1, quiet = False):
    _validate_inputs(fl, fh, win_length_sec, overlap, upsample_ratio)

    # arguments to be provided to xcorr_one_day.
    # 'quiet' is provided to loop_through_days directly, as well as to xcorr_one_day via args.
    args = {'fl':fl, 'fh':fh, 'window_length_sec':win_length_sec, 'overlap':overlap, 'upsample_ratio':upsample_ratio, 'quiet':quiet}
    
    ## separate validations are needed here (to prevent awkward-to-resolve errors in loop_through_days)

    return loop_through_days(xcorr_one_day, files, t1, t2, IDs, exclude_IDs = exclude_IDs, quiet = quiet, args = args)
                             
#########################################################
#########################################################
def invert_for_slowness(xcorr_df, locations):
    """
    locations: pd.DataFrame with columns x, y, and network/station/location (output of get_coordinates)
    """
    data_short_IDs = []
    for key in xcorr_df.keys():
        if re.search('rms', key):
            data_short_IDs.append(key.split('_')[1])
            
    locations['ID'] = [f'{locations.network[i]}.{locations.station[i]}.{locations.location[i]}' for i in range(locations.shape[0])]
    keep_indices = np.where(np.isin(locations.ID, data_short_IDs))[0]
    locations = locations.iloc[keep_indices, :]
    locations = locations.sort_values('ID', ignore_index = True)

    ## solve linear system G . s = t: G is x/y distances, s is slowness, and t is observed time lags
    G = []
    lag_keys = []
    #for key in xcorr_df.keys():
    #    if re.search('lag', key):
    for i, short_ID_1 in enumerate(data_short_IDs):
        short_ID_2 = data_short_IDs[(i+1) % len(data_short_IDs)]
        lag_keys.append(f'lag_{short_ID_1}_{short_ID_2}')
        #(short_ID_1, short_ID_2) = key.split('_')[1:]
        index1 = np.where(locations.ID == short_ID_1)[0][0]
        index2 = np.where(locations.ID == short_ID_2)[0][0]
        G.append([locations.loc[index1, 'x'] - locations.loc[index2, 'x'],
                  locations.loc[index1, 'y'] - locations.loc[index2, 'y']])

    G = np.array(G)[:-1,:] # drop the last row because it's linearly dependent on the others
    
    ## Use the generalized inverse in case there are more than 3 sensors: (GTG)^-1 . GT . t = H . t = s
    ## @ is matrix multiplication symbol
    H = np.linalg.inv(G.T @ G) @ G.T

    ## loop through all the time windows and make a slowness vector for each
    n_windows = xcorr_df.shape[0]
    backazimuth = np.zeros(n_windows)
    slowness = np.zeros(n_windows)
    for i in range(n_windows):
        lags = np.array([xcorr_df[lag_key][i] for lag_key in lag_keys[:-1]])
        slowness_vector = H @ lags
        backazimuth[i] = np.arctan2(-slowness_vector[0], -slowness_vector[1]) * 180/np.pi
        slowness[i] = np.sqrt(np.sum(slowness_vector**2))
    xcorr_df['backazimuth'] = backazimuth
    xcorr_df['slowness'] = slowness
    return xcorr_df
#########################################################
#########################################################
    

def loop_through_days(function, filenames, t1 = '1970-01-01', t2 = '9999-12-31', IDs = None, exclude_IDs = None, quiet = False, args = {}):
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
    if len(filenames) == 0:
        raise Exception('No files found; check data files')
    for filename in filenames:
        try:
            st = obspy.read(filename, headonly=True)
        except:
            print(f'skipping unreadable file {filename}')
            continue
        for tr in st:
            file_metadata['filename'].append(filename)
            file_metadata['t1'].append(tr.stats.starttime)
            file_metadata['t2'].append(tr.stats.endtime)
            file_metadata['ID'].append(tr.id)
    
    file_metadata_df = pd.DataFrame.from_dict(file_metadata)

    t1 = np.max([t1, file_metadata_df.t1.min()])
    t2 = np.min([t2, file_metadata_df.t2.max()])
    IDs = sorted(_check_input_IDs(file_metadata_df, IDs, exclude_IDs))

    print('Processing data from channels ' + ', '.join(IDs))

    rows_to_keep = np.where((file_metadata_df.t2 >= t1) & \
                    (file_metadata_df.t1 <= t2) & \
                    np.isin(file_metadata_df.ID, IDs))[0]
    file_metadata_df = file_metadata_df.iloc[rows_to_keep,:]
    file_metadata_df.sort_values('t1', ignore_index=True, inplace=True)

    if file_metadata_df.shape[0] == 0:
        raise Exception('No data files fit t1/t2/IDs criteria; check those inputs')
    
    ## loop through days. be careful to avoid funny business with leap seconds.
    day_starts = [t1]
    while True:
        new_time = (day_starts[-1]+86400+2).replace(hour=0, minute=0, second=0) # add 24 hours + 2 sec, then round down
        if new_time < t2:
            day_starts.append(new_time)
        else:
            break
    day_starts = np.array(day_starts)

    day_ends = [day_starts[0]+86400]
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
        if not quiet:
            print(f'Processing from {day_start} to {day_end}')
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
        #for i, tr in enumerate(st):
        #    if tr.id not in IDs:
        #        st.pop(i) 
        st = obspy.Stream([tr for tr in st if tr.id in IDs])
        
        ## finally, apply whatever function you have to the data.
        ## 'function' must accept two inputs: an obspy.Stream with data,
        ## and a dictionary with function-specific arguments.
        ## the function will return a pd.DataFrame. append to a list, then merge at the end
        try:
            day_output = function(st.slice(day_start, day_end), args)
            output_list.append(day_output)
        except:
            print(f'Error on day {t1.isoformat()}, skipping')
            print(traceback.format_exc())
    ## done looping. merge the output and return.
    if len(output_list) == 0:
        return None
    else:
        return pd.concat(output_list, ignore_index = True)
        
def _validate_inputs(fl, fh, win_length_sec, overlap, upsample_ratio):
    try:
        if (fl < 0) or (fh < 0) or (fl >= fh):
            raise Exception(f'fl ({fl}) and fh ({fh}) must both be non-negative numbers or NaN, and fh > fl')
    except:
        raise Exception(f'fl ({fl}) and fh ({fh}) must both be non-negative numbers or NaN, and fh > fl')

    try:
        if win_length_sec <= 0:
            raise Exception(f'win_length_sec ({win_length_sec}) must be a positive number')
    except:
        raise Exception(f'win_length_sec ({win_length_sec}) must be a positive number')

    try:
        if not ((overlap >= 0) and (overlap < 1)):
            raise Exception('overlap must be a non-negative number strictly less than 1')
    except:
        raise Exception('overlap must be a non-negative number strictly less than 1')

    try:
        if not (upsample_ratio >= 1):
            raise Exception('upsample_ratio must be at least 1')
    except:
        raise Exception('upsample_ratio must be at least 1')

    return
    
def xcorr_one_day(st, args = {}):
    ## check inputs and implement defaults
    fl = args.get('fl') if args.get('fl') is not None else 5
    fh = args.get('fh') if args.get('fh') is not None else 40
    win_length_sec = args.get('win_length_sec') if args.get('win_length_sec') is not None else 5
    overlap = args.get('overlap') if args.get('overlap') is not None else 0
    upsample_ratio = args.get('upsample_ratio') if args.get('overlap') is not None else 1
    _validate_inputs(fl, fh, win_length_sec, overlap, upsample_ratio)

    if upsample_ratio != 1:
        st = upsample_stream(st, args['upsample_ratio'])
        
    st.detrend('linear')

    ## de-step the beginning of the trace to prevent filter artifacts
    for tr in st:
        tr.data -= tr.data[0]
    st.filter('bandpass', freqmin = fl, freqmax = fh)
    output = apply_function_windows(st, xcorr_function, win_length_sec, overlap, args)

    ## reformat UTCDateTimes as string
    #output['t_mid'] = [t.isoformat() for t in output['t_mid']]
    return pd.DataFrame.from_dict(output)

def xcorr_function(st, args):
    maxshift_seconds = args.get('maxshift_seconds')
    if maxshift_seconds is None:
        maxshift_seconds = 1

    quiet = args.get('quiet')
    if quiet is None:
        quiet = False
    #if not quiet:
    #    print(st)
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
        ID1 = f'{tr1.stats.network}.{tr1.stats.station}.{tr1.stats.location}'
        ID2 = f'{tr2.stats.network}.{tr2.stats.station}.{tr2.stats.location}'
        output_dict[f'rms_{ID1}'] = tr1.std()
        pair_name = f'{ID1}_{ID2}'
        shift, value = xcorr_max(correlate(tr1.data, tr2.data, int(np.round(maxshift_seconds / dt))), abs_max = False)
        output_dict[f'lag_{pair_name}'] = shift * dt
        output_dict[f'r_{pair_name}'] = value
        consistency += shift/len(st)
        output_dict['mean_coef'] += value/len(st)
    ## consistency: allow up to one sample error per cross-correlation
    output_dict['consistency'] = np.abs(consistency) <= len(st)
    return output_dict



def apply_function_windows(st, f, win_length_sec, overlap = 0.5, args = {}):
    """
    Run an analysis (or suite of analyses) on overlapping windows for some dataset
    
    Parameters:
    -----------
    st : obspy.Stream
    Stream including data to divide into windows and analyze

    f : function
    Accepts variable "st" (obspy.Stream) and "args" (dict), returns dictionary of results

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
    for i in range(num_windows):
        win_start = t1 + i*(data_length_sec - win_length_sec) / (num_windows-1)
        st_tmp = st.slice(win_start-eps, win_start + win_length_sec - eps, nearest_sample = False)
        win_dict = f(st_tmp, args)
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


def _check_input_IDs(file_metadata_df, IDs, exclude_IDs):
    """
    Validate station IDs provided by user
    """
    found_IDs = _unique(file_metadata_df.ID)
    ## if no IDs are provided by the user, assume that all IDs in the data are allowed
    if (IDs is None) or (len(IDs) == 0):
        output_IDs = found_IDs
    else:
        ## if the user does provide IDs, we need to validate all of them
        output_IDs = []
        for ID in IDs:
            ID = str(ID)
            ## first, check to see if the ID is formatted exactly per naming convention
            ## if yes, pass it unmodified.
            ## https://ds.iris.edu/ds/nodes/dmc/data/formats/seed/
            if re.match(r'\w{0,2}\.\w{1,5}\.\w{0,2}\.\w{3}', ID):
                output_IDs.append(ID)

            ## Next, search for IDs in the data that match the provided ID.
            ## There's no reason the user can't provide a regex here.
            else:
                for found_ID in found_IDs:
                    if re.search(ID, found_ID):
                        output_IDs.append(found_ID)
        output_IDs = _unique(output_IDs)
    if exclude_IDs is not None:
        remaining_IDs = []
        for i, ID in enumerate(output_IDs):
            exclude_this_ID = False
            for exclude_ID in exclude_IDs:
                if re.search(exclude_ID, ID):
                    exclude_this_ID = True
            if not exclude_this_ID:
                remaining_IDs.append(ID)
        output_IDs = remaining_IDs
    return(output_IDs)

def get_coordinates(x, y = None):
    """
    Finds sensor coordinates from various inputs
    
    Parameters:
    -----------
    x: either an array of x coordinates, obspy.Stream with trace.stats['coords'], or obspy.Inventory
    y: either an array of y coordinates (if x is an array of x coordinates), or None

    Returns:
    --------
    pandas.DataFrame with x, y, z, network, station, location fields. x and y are in km, z is in m.
    """
    import obspy.signal.array_analysis
    if type(x) is obspy.Stream:
        try:
            ## Stream coordinates can either be lon/lat/z or x/y/z. x and y are km, z is m.
            ## This line will work if x has lon/lat/z coordinates, and will raise an exception
            ## if x has x/y/z coordinates.
            geometry = obspy.signal.array_analysis.get_geometry(x, coordsys = 'lonlat',
                                                              return_center = True)
            coords = {'x':geometry[:-1,0],
                      'y':geometry[:-1,1],
                      'z':geometry[:-1,2] + geometry[-1,2]} # last row in 'geometry' is coordinates of reference point
            # https://docs.obspy.org/packages/autogen/obspy.signal.array_analysis.get_geometry.html
        except:
            ## If we're here, then x is a stream with x/y/z coordinates. Extract them directly.
            ## We don't want to use get_geometry because it will pick a new center and shift the
            ## coordinates accordingly.
            coords = {'x': np.array([tr.stats.coordinates['x'] for tr in x]),
                      'y': np.array([tr.stats.coordinates['y'] for tr in x]),
                      'z': np.array([tr.stats.coordinates['elevation'] for tr in x])}

        coords['network'] = [tr.stats.network for tr in x]
        coords['station'] = [tr.stats.station for tr in x]
        coords['location'] = [tr.stats.location for tr in x]
    elif type(x) is obspy.Inventory:
        contents = x.get_contents()['channels']
        lats = [x.get_coordinates(s)['latitude'] for s in contents]
        lons = [x.get_coordinates(s)['longitude'] for s in contents]
        zz = [x.get_coordinates(s)['elevation'] for s in contents]
        xx = np.zeros(len(lats))
        yy = np.zeros(len(lats))
        for i, (lat, lon) in enumerate(zip(lats, lons)):
            xx[i], yy[i] = obspy.signal.util.util_geo_km(np.mean(lons), np.mean(lats), lon, lat)
        coords = {'x': xx, 'y': yy, 'z': zz,
                  'network': [string.split('.')[0] for string in contents],
                  'station': [string.split('.')[1] for string in contents],
                  'location': [string.split('.')[2] for string in contents]}
    else:
        coords = {'x':x,
                  'y':y,
                  'z':None,
                  'network':None,
                  'station':None,
                  'location':None}
    return pd.DataFrame(coords)


def upsample_stream(st, N):
    for tr in st:
        t_in = np.arange(tr.stats.npts)
        t_out = np.arange(N * (tr.stats.npts - 1)) / N
        spline = scipy.interpolate.CubicSpline(t_in, tr.data)
        tr.data = spline(t_out)
        tr.stats.delta /= N
    return st

