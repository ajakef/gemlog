import pdb
import warnings
import numpy as np
from numpy import NaN, Inf
import os, glob, csv, time, scipy
import pandas as pd
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
import obspy
#import datetime
import sys

_debug = False

class EmptyRawFile(Exception):
    """Raised when the input file is empty"""
    pass
class CorruptRawFile(Exception):
    """Raised when the input file cannot be read for some reason"""
    pass

class CorruptRawFileNoGPS(CorruptRawFile):
    """Raised when the input file does not contain any GPS data"""
    pass

class MissingRawFiles(Exception):
    """Raised when no input files are readable"""
    pass


def _breakpoint():
    if _debug: # skip if we aren't in debug mode
        pdb.set_trace()
#####################
def convert(rawpath = '.', convertedpath = 'converted', metadatapath = 'metadata', \
            metadatafile = '', gpspath = 'gps', gpsfile = '', t1 = -Inf, t2 = Inf, nums = NaN, \
            SN = '', bitweight = NaN, units = 'Pa', time_adjustment = 0, blockdays = 1, \
            file_length_hour = 1, station = '', network = '', location = '', output_format = 'MSEED'):
    """
    Read raw Gem files, interpolate them, and write output files in miniSEED or SAC format.

    Parameters
    ----------
    rawpath : str, default '.' 
        Path of folder containing raw Gem data files to read.

    convertedpath : str, 'converted'
        Path of folder where converted data files should be written. If this folder does 
        not exist, convert will try to create it.
    
    metadatapath : str, 'metadata'
        Path of folder where metadata files should be written. If this folder does 
        not exist, convert will try to create it.
    
    metadatafile : str
        File name for metadata file. Default is of the form 'XXXmetadata_NNN.txt', where
        'XXX' is the Gem's serial number and 'NNN' is the file index (000 at first, incrementing by one for each subsequent calculation).

    gpspath : str, 'gps'
        Path of folder where gps files should be written. If this folder does 
        not exist, convert will try to create it.
    
    gpsfile : str
        File name for gps file. Default is of the form 'XXXgps_NNN.txt', where
        'XXX' is the Gem's serial number and 'NNN' is the file index (000 at first, incrementing by one for each subsequent calculation).
    
    t1 : float or obspy.UTCDateTime, default -np.inf
        Start time for the conversion. If float, the number of seconds 
        since the epoch (1970-01-01 00:00:00 UTC). 

    t2 : float or obspy.UTCDateTime, default np.inf
        Stop time for the conversion. If float, the number of seconds 
        since the epoch (1970-01-01 00:00:00 UTC). 

    nums : list or np.array of integers
        Numbers of raw Gem files to read. By default, it reads all files
        in rawpath for the specified serial number.
    
    SN : str
        One Gem serial number to read and convert. Use a loop to convert
        multiple Gems.

    bitweight : float
        The value of each count when converting between counts and other 
        units (typically Pascals, possibly Volts). By default, it looks
        up the correct bitweight given the Gem's serial number and gain 
        configuration. Leave this blank unless the Gem has been modified
        in a way that changes the bitweight.

    units : str, default 'Pa'
        Desired output units. Options are 'Pa' (Pascals), 'V' (Volts), 
        or 'counts'.

    time_adjustment : float, default 0
        Amount to shift the sample times by in case of timing error, 
        possibly due to leap second issues.

    blockdays : float, default 1
        Number of days worth of data to read in and convert at a time.
        If you have memory issues when converting data, try setting this
        to a smaller value.

    file_length_hour : float, default 1
        Length of the output files in hours.
            
    station : str
        Name of the station (up to five characters) to assign to the 
        data. If not provided, uses the Gem's serial number as the 
        station ID.

    network : str
        Two-character name of the sensor network. Leaving this blank is
        normally fine in subsequent data processing.

    location : str
        Two-character location code for this Gem. Leaving this blank is 
        usually fine in subsequent data processing.
    
    output_format : str, default 'MSEED'
        Output file format. Currently, formats 'MSEED' and 'SAC' are 
        supported.

    Returns
    -------
    None, writes output files only (converted, metadata, and gps)

    Note
    ----
    All sample times involving the Gem (and most other passive 
    seismic/acoustic data) are in UTC; time zones are not supported.
    """
    file_length_sec = 3600 * float(file_length_hour)
    ## bitweight: leave blank to use default (considering Gem version, config, and units). This is preferred when using a standard Gem (R_g = 470 ohms)
    
    ## make sure the raw directory exists and has real data
    if not os.path.isdir(rawpath):
        raise Exception('Raw directory ' + rawpath + ' does not exist')
    if len(glob.glob(rawpath + '/FILE' +'[0-9]'*4 + '.???')) == 0:
        raise Exception('No data files found in directory ' + rawpath)

    ## make sure bitweight is a scalar
    if((type(nums) == type(1)) or (type(nums) == type(1.0))):
        nums = np.array([nums])
    else:
        nums = np.array(nums)

    ## if 'nums' is default, convert all the files in this directory
    if((len(nums) == 0) or np.isnan(nums[0])):
        if(True or len(SN) == 1): # trying to check that SN is a scalar
            fn = glob.glob(rawpath + '/FILE' +'[0-9]'*4 + '.' + SN) + \
                 glob.glob(rawpath + '/FILE' +'[0-9]'*4 + '.TXT')
        else:
            fn = glob.glob(rawpath + '/FILE' +'[0-9]'*4 + '.???')
        nums = np.array([int(x[-8:-4]) for x in fn]) # "list comprehension"

    ## Catch if rawpath doesn't contain any files from SN. This won't catch files ending in TXT.
    if len(nums) == 0:
        raise Exception('No data files for SN "' + SN + '" found in raw directory ' + rawpath)
    nums.sort()

    ## start at the first file in 'nums'
    n1 = np.min(nums)
  
    ## read the first set of up to (24*blockdays) files
    L = _new_gem_var()
    while((L['data'].count() == 0) & (n1 <= max(nums))): ## read sets of files until we get one that isn't empty
        nums_block = nums[(nums >= n1) & (nums < (n1 + (12*blockdays)))] # files are 2 hours, so 12 files is 24 hours
        L = read_gem(nums_block, rawpath, SN = SN, network = network, station = station, location = location)
        n1 = n1 + (12*blockdays) # increment file number counter

    p = L['data']
    
    ## if bitweight isn't set, use the default bitweight for the logger version, config, and units
    if(np.isnan(bitweight)):
        if(units == 'Pa'):
            bitweight = L['header']['bitweight_Pa'][0]
        elif(units == 'V'):
            bitweight = L['header']['bitweight_V'][0]
        elif (units == 'Counts') | (units == 'counts'):
            bitweight = 1
      
    ## if not specified, define t1 as the earliest integer-second time available
    if(np.isinf(float(t1))):
        p = L['data']
        p.merge()
        t1 = p[0].stats.starttime
        t1 = obspy.core.UTCDateTime(np.ceil(float(t1)))

    if(np.isinf(float(t2))):
        t2 = obspy.core.UTCDateTime.strptime('9999-12-31 23:59:59', '%Y-%m-%d %H:%M:%S') # timekeeping apocalypse
  
    wsn = 0
    while(len(SN) < 3): # take the first non-NA SN. This is important because there can be blank files in there.
        wsn = wsn+1
        SN = L['header']['SN'][wsn]
  
    ## set up the gps and metadata files. create directories if necessary
    if(len(gpsfile) == 0):
        if(not os.path.isdir(gpspath)):
            try:
                os.makedirs(gpspath) # makedirs vs mkdir means if gpspath = dir1/dir2, and dir1 doesn't exist, that dir1 will be created and then dir1/dir2 will be created
            except:
                print('Failed to make directory ' + gpspath)
                sys.exit(2)
        gpsfile = _make_filename(gpspath, SN, 'gps')

  
    if(len(metadatafile) == 0):
        if(not os.path.isdir(metadatapath)):
            try:
                os.makedirs(metadatapath)
            except:
                print('Failed to make directory ' + metadatapath)
                sys.exit(2)
        metadatafile = _make_filename(metadatapath, SN, 'metadata')
  
    ## if the converted directory does not exist, make it
    if(not os.path.isdir(convertedpath)):
        try:
            os.makedirs(convertedpath)
        except:
            print('Failed to make directory ' + convertedpath)
            sys.exit(2)
  
    ## start metadata and gps files
    metadata = L['metadata']   
    gps = L['gps']
    metadata.to_csv(metadatafile, index=False) ## change to metadata format. need to make ScanMnetadata compatible with both

    wgps = (gps['t'] > (t1 - 1)) 
    if(len(wgps) > 0):
        gps[wgps].to_csv(gpsfile, index=False)

    hour_to_write = max(t1, p[0].stats.starttime)
    hour_to_write = _write_hourlong_mseed(p, hour_to_write, file_length_sec, bitweight, convertedpath, output_format=output_format)
    
    ## read sets of (12*blockdays) files until all the files are converted
    while(True):
        ## check to see if we're done
        if(n1 > np.max(nums)):# & len(p) == 0):
            break # out of raw data to convert
        if((t1 > t2) & (not np.isnan(t1 > t2))):
            break # already converted the requested data
        ## load new data if necessary
        tt2 = min(t2, _trunc_UTCDateTime(t1, 86400*blockdays) + 86400*blockdays)
        while((p[-1].stats.endtime < tt2) & (n1 <= max(nums))):
            L = read_gem(nums[(nums >= n1) & (nums < (n1 + (12*blockdays)))], rawpath, SN = SN)
            n1 = n1 + (12*blockdays) # increment file counter

            if(len(L['data']) == 0):
                continue # skip ahead if there aren't any readable data files here

            ## process newly-read data
            if(any(L['header'].SN != SN) | any(L['header'].SN.apply(len) == 0)):
                #_breakpoint()
                w = (L['header'].SN != SN) | (L['header'].SN.apply(len) == 0)
                #print('Wrong or missing serial number(s): ' + L['header'].SN[w] + ' : numbers ' + str(nums[np.logical_and(nums >= n1, nums < (n1 + (12*blockdays)))][w]))
                for i in w:
                    print('Wrong or missing serial number: ' + L['header'].file[i])

            #pdb.set_trace()
            if(len(L['data']) > 0):
                p = p + L['data']
                p.merge()
            #print(p)
                
            ## start metadata and gps files
            metadata = L['metadata']   
            gps = L['gps']
            metadata.to_csv(metadatafile, index=False, mode='a', header=False)
            wgps = (gps['t'] > (t1 - 1)) 
            ## update the gps file
            if(len(wgps) > 0):
                gps.to_csv(gpsfile, index=False, mode='a', header=False)
                
        ## run the conversion and write new converted files
        while((hour_to_write + file_length_sec) <= p[-1].stats.endtime):
            hour_to_write = _write_hourlong_mseed(p, hour_to_write, file_length_sec, bitweight, convertedpath, output_format=output_format)
            
        ## update start time to convert
        p.trim(hour_to_write, t2)
        t1 = _trunc_UTCDateTime(tt2+(86400*blockdays) + 1, 86400*blockdays)
    ## done reading new files. write what's left and end.
    while((hour_to_write <= p[-1].stats.endtime) & (len(p) > 0)):
        hour_to_write = _write_hourlong_mseed(p, hour_to_write, file_length_sec, bitweight, convertedpath, output_format=output_format)
        p.trim(hour_to_write, t2)
        p = p.split()
        if(len(p) > 0):
            hour_to_write = p[0].stats.starttime
        else:
            break

Convert = convert # alias; v1.0.0

def _write_hourlong_mseed(p, hour_to_write, file_length_sec, bitweight, convertedpath, hour_end = np.nan, output_format='mseed'):
    #pdb.set_trace()
    if(np.isnan(hour_end)):
        hour_end = _trunc_UTCDateTime(hour_to_write, file_length_sec) + file_length_sec
    pp = p.copy()
    pp.trim(hour_to_write, hour_end)
    pp = pp.split() ## in case of data gaps ("masked arrays", which fail to write)
    #_breakpoint()
    for tr in pp:
        tr.stats.calib = bitweight
        fn = _make_filename_converted(tr, output_format)
        if(len(tr) > 0):
            print(tr)
            if(output_format.lower() == 'wav'):
                ## this is supposed to work for uint8, int16, and int32, but actually only works for uint8. obspy bug?
                    tr.data = np.array(tr.data, dtype = 'uint8')# - np.min(tr[i].data)
                    tr.write(convertedpath +'/'+ fn, format = 'WAV', framerate=7000, width=1) 
            else:
                tr.write(convertedpath +'/'+ fn, format = output_format, encoding=10) # encoding 10 is Steim 1
        #mseed_core._write_mseed(pp, convertedpath +'/'+ fn, format = 'MSEED', encoding=10)

    hour_to_write = hour_end
    return hour_to_write
    
####################################

def _trunc_UTCDateTime(x, n=86400):
    return obspy.core.UTCDateTime(int(float(x)/n)*n)#, origin='1970-01-01')

def _make_filename(dir, SN, dirtype):
    n = 0
    fn = dir + '/' + SN + dirtype + '_'+ f'{n:03}' + '.txt'
    while(os.path.exists(fn)):
        n = n + 1
        fn = dir + '/' + SN + dirtype + '_' + f'{n:03}' + '.txt'

    return fn


def _make_filename_converted(pp, output_format):
    t0 = pp.stats.starttime
    ## colons separating H:M:S would be more readable, but are not allowed in Windows filenames
    return f'{t0.year:04}' + '-' +f'{t0.month:02}' + '-' +f'{t0.day:02}' + 'T' + f'{t0.hour:02}' + '_' + f'{t0.minute:02}' + '_' + f'{t0.second:02}' + '.' + pp.id + '.' + output_format.lower()
#import pdb

def _new_gem_var():
    tr = obspy.Trace()
    tr.stats.delta = 0.01
    gps = pd.DataFrame(columns=['year', 'date', 'lat', 'lon'])
    metadata = pd.DataFrame(columns=['millis', 'batt', 'temp', 'A2', 'A3', \
                                     'maxWriteTime', 'minFifoFree', 'maxFifoUsed', \
                                     'maxOverruns', 'gpsOnFlag', 'unusedStack1',\
                                     'unusedStackIdle', 't'])
    output = {'data': tr,
              'metadata': metadata,
              'gps': gps
    }
    return output


def make_db(path, pattern = '*', savefile = './DB.csv'):
    """Create a database summarizing a set of converted data files.

    Parameters
    ----------
    path : str
        Path to folder containing converted data files to summarize.
    
    pattern : str, default '*'
        Glob-type pattern for selecting converted data files to summarize

    savefile : str, default './DB.csv'
        File name where database should be written.

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
        if((count % 100) == 0):
            print(str(count) + ' of ' + str(len(files)))
        count = count + 1
    DB = pd.concat(DB)
    DB.to_csv(savefile)
    return(DB)

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


################################################
def _read_SN(fn):
    SN_line = pd.read_csv(fn, delimiter = ',', skiprows = 4, nrows=1, dtype = 'str', names=['s', 'SN'])
    SN = SN_line['SN'][0]
    return SN

def _read_format_version(fn):
    versionLine = pd.read_csv(fn, delimiter = ',', nrows=1, dtype = 'str', names=['s'])
    version = versionLine['s'][0][7:]
    return version
    
def _read_config(fn):
    config = pd.Series({'gps_mode': 1,
                    'gps_cycle' : 15,
                    'gps_quota' : 20,
                    'adc_range' : 0,
                    'led_shutoff' : 0,
                    'serial_output' : 0}) ## default config: it's fairly safe to use this as the default because any other configuration would require ...
    for j in range(10):
        line = pd.read_csv(fn, skiprows = j+1, nrows=1, delimiter = ',', dtype = 'str', names = ['na', 'gps_mode','gps_cycle','gps_quota','adc_range','led_shutoff','serial_output'])
        if line.iloc[0,0] == 'C':
            #config = line.iloc[0,1:]
            config = {key:int(line[key]) for key in list(line.keys())[1:]}
            break
    return config

def _fn2nums(fn_list):
    nums = []
    for i, fn in enumerate(fn_list):
        nums[i] = int(fn[-8:-5])
    return nums


def _find_nonmissing_files(path, SN, nums):
    ## list all Gem files in the path
    fnList = glob.glob(path + '/' + 'FILE[0-9][0-9][0-9][0-9].[0-9][0-9][0-9]')
    fnList.sort()
    fnList = np.array(fnList)
    
    ## find out what all the files' SNs are
    ext = np.array([x[-3:] for x in fnList])
    for i in range(len(ext)):
        if ext[i] == 'TXT':
            ext[i] = _read_SN(fnList[i])
    ## check the files for SN and num
    goodFnList = []
    for i, fn in enumerate(fnList):
        fnNum = int(fn[-8:-4])
        fnSN = ext[i]
        if (fnNum in nums) & (fnSN == SN):
            goodFnList.append(fn)
    if len(goodFnList) == 0:
        print('No good data files found for specified nums and SN ' + SN)
        return []
        ## fix this to be an exception or warning?
    ## make sure they aren't empty
    goodNonemptyFnList = []
    for fn in goodFnList:
        if os.path.getsize(fn) > 10: # to be safe, anything under 10 bytes is treated as empty 
            goodNonemptyFnList.append(fn)
    if(len(goodNonemptyFnList) == 0):
        ## warning
        print('No non-empty files')
        return []
    if(len(goodNonemptyFnList) < len(goodFnList)):
        print('Some files are empty, skipping')
    return goodNonemptyFnList    



def _unwrap_millis(new, old, maxNegative = 2**12, rollover = 2**13):
    ## maxNegative is the greatest allowable negative difference.
    ## negative differences can happen between data and GPS lines, or between metadata and data lines.
    return old + ((new - (old % rollover) + maxNegative) % rollover) - maxNegative

def _check_gps(line): # return True if GPS line is good
    #G,msPPS,msLag,yr,mo,day,hr,min,sec,lat,lon
    return not ((line[8] == 0) or (line[8] > 90) or (line[8] < -90) or # lat
                (line[9] == 0) or (line[9] > 180) or (line[9] < -180) or # lon
                (line[1] > 1000) or (line[1] < 0) or # lag
                (line[2] > 2040) or (line[2] < 2014) or # year
                (line[3] > 12) or (line[3] < 1) or # month
                (line[4] > 31) or (line[4] < 1) or # day
                (line[5] > 24) or (line[5] < 0) or # hour
                (line[6] > 60) or (line[6] < 0) or # minute
                (line[7]>60) or (line[7]<0) or (line[7]!=np.round(line[7]))) # second


def _make_gps_time(line):
    try:
        return obspy.UTCDateTime(int(line[2]), int(line[3]), int(line[4]), int(line[5]), int(line[6]), int(line[7]))
    except:
        return np.NaN

def _millis_to_time(L):
    G = L['gps']
    D = L['data']
    coefficients = np.polyfit(G.msPPS, G.t, 3)    
    #print(coefficients)
    pf = np.poly1d(coefficients)
    return pf


def _read_with_cython(filename, offset=0, require_gps = True):
    """
    Read a Gem logfile.

    Parameters
    ----------
    filename : str or pathlib.Path
        Filepath of a file containing data to read.
    offset : int, default 0
        A timing offset to include in the millisecond timestamp values.

    Returns
    -------
    dict of dataframes
        - data: the analog readings and associated timings
        - metadata: datalogger metadata
        - gps: GPS timing and location values
    """
    try:
        from gemlog.parsers import parse_gemfile
    except ImportError:
        raise ImportError(
            "gemlog's C-extensions are not available. Reinstall gemlog with "
            "C-extensions to use this function."
        )

    # use cythonized reader file instead of pd.read_csv and slow string ops
    values, types, millis = parse_gemfile(str(filename).encode('utf-8'))
    if values.shape[0] == 0:
        raise EmptyRawFile(filename)
    if (b'G' not in types) and require_gps:
        raise CorruptRawFileNoGPS(filename)
    df = pd.DataFrame(values, columns=range(2, 13))
    # note that linetype has type bytes here, not str like in the pandas func
    df['linetype'] = types
    df['millis-sawtooth'] = millis
    return _process_gemlog_data(df, offset)


def _read_with_pandas(filename, offset=0, require_gps = True):
    # skiprows is important so that the header doesn't force dtype=='object'
    # the C engine for pd.read_csv is fast but crashes sometimes. Use the python engine as a backup.
    try:
        df = pd.read_csv(filename, names=range(13), low_memory=False, skiprows=6)
    except Exception:
        try:
            df = pd.read_csv(filename, names=range(13), engine='python', skiprows=6,
                             error_bad_lines = False, warn_bad_lines = False)
        except:
            raise CorruptRawFile(filename)
    if df.shape[0] == 0:
        raise EmptyRawFile(filename)

    try:
        df['linetype'] = [value[0] for value in df[0]] # exception if any 'value' is not a string
    except:
        raise CorruptRawFile(filename)
    if ('G' not in set(df.loc[:,'linetype'])) and require_gps:
        raise CorruptRawFileNoGPS(filename)
    
    df = df.loc[df.loc[:,'linetype'].isin(['D', 'M', 'G']), :]
    #df = df.loc[df['linetype'].isin(['D', 'M', 'G']), :]
    ## most of the runtime is before here
    # unroll the ms rollover sawtooth
    df['millis-sawtooth'] = np.where(df['linetype'] == 'D',df[0].str[1:],df[1]).astype(int)
    return _process_gemlog_data(df, offset)


def _read_single_v0_9(filename, offset=0, require_gps = True):
    """
    Read a Gem logfile.

    Parameters
    ----------
    filename : str or pathlib.Path
        Filepath of a file containing data to read.

    offset : int, default 0
        A timing offset to include in the millisecond timestamp values.

    Returns
    -------
    dict of dataframes

        - data: the analog readings and associated timings
        - metadata: datalogger metadata
        - gps: GPS timing and location values
    """
    # Try each of the three file readers in order of decreasing speed but
    # probably increasing likelihood of success.

    readers = [
        _read_with_cython, _read_with_pandas, _slow__read_single_v0_9
    ]
    for reader in readers:
        try:
            output = reader(filename, offset, require_gps)
            if (len(output['gps'].lat) == 0) and require_gps:
                raise CorruptRawFileNoGPS(filename)
            return output
        except (EmptyRawFile, FileNotFoundError, CorruptRawFileNoGPS):
            # If the file is definitely not going to work, exit early and
            # re-raise the exception that caused the problem
            raise
        except Exception:
            pass

    raise CorruptRawFile(filename)


def _process_gemlog_data(df, offset):

    # unroll the ms rollover sawtooth
    df['millis-stairstep'] = (df['millis-sawtooth'].diff() < -2**12).cumsum()
    df['millis-stairstep'] -= (df['millis-sawtooth'].diff() > 2**12).cumsum()
    df['millis-stairstep'] *= 2**13
    df['millis-corrected'] = df['millis-stairstep'] + df['millis-sawtooth']
    first_millis = df['millis-corrected'].iloc[0]
    df['millis-corrected'] += (
        (offset-first_millis)
        + ((first_millis-(offset % 2**13)+2**12) % 2**13)
        - 2**12
    )
    # groupby is somewhat faster than repeated subsetting like
    # df.loc[df['linetype'] == 'D', :], ...
    grouper = df.groupby('linetype')
    # the python-based reader uses strings for linetype but the cython version
    # uses bytes.  figure out which one we need:
    if isinstance(df['linetype'].iloc[0], bytes):
        Dkey = b'D'
        Gkey = b'G'
        Mkey = b'M'
        data_col = 2
    else:
        Dkey, Gkey, Mkey = 'DGM'
        data_col = 1
    D = grouper.get_group(Dkey)
    G = grouper.get_group(Gkey)
    M = grouper.get_group(Mkey)

    # pick out columns of interest and rename
    D_cols = ['msSamp', 'ADC']
    G_cols = ['msPPS', 'msLag', 'year', 'month', 'day', 'hour', 'minute',
              'second', 'lat', 'lon']
    M_cols = ['millis', 'batt', 'temp', 'A2', 'A3',
              'maxWriteTime', 'minFifoFree', 'maxFifoUsed',
              'maxOverruns', 'gpsOnFlag', 'unusedStack1', 'unusedStackIdle']

    # column names are currently integers (except for the calculated cols)
    D = D[['millis-corrected', data_col]]
    D.columns = D_cols
    G = G[['millis-corrected'] + list(range(2, len(G_cols)+1))]
    G.columns = G_cols
    M = M[['millis-corrected'] + list(range(2, len(M_cols)+1))]
    M.columns = M_cols

    # now that columns aren't mixed dtype anymore,
    # convert to numeric where possible
    D = D.apply(pd.to_numeric)
    G = G.apply(pd.to_numeric)
    M = M.apply(pd.to_numeric)

    # filter bad GPS data and combine into datetimes
    valid_gps = _valid_gps(G)
    G = G.loc[valid_gps, :]

    def make_gps_time(row):
        try:
            return obspy.UTCDateTime(*row)
        except Exception:
            return np.NaN

    G['t'] = G.iloc[:, 2:8].astype(int).apply(make_gps_time, axis=1)

    # process data (version-dependent)
    D['ADC'] = D['ADC'].astype(float).cumsum()

    return {'data': np.array(D),
            'metadata': M.reset_index().astype('float'),
            'gps': G.reset_index().astype('float')}


def _valid_gps(G):
    # vectorized GPS data validation
    # basic lower and upper bounds:
    limits = {
        'lat': (-90, 90),
        'lon': (-180, 180),
        'msLag': (0, 1000),
        'year': (2014, 2040),
        'month': (1, 12),
        'day': (1, 31),
        'hour': (0, 24),
        'minute': (0, 60),
        'second': (0, 60),
    }
    bad_gps = False
    for key, (lo, hi) in limits.items():
        data = G[key]
        bad_gps |= (data < lo) | (data > hi)

    # custom limits:
    bad_gps |= (
        (G['lat'] == 0) |
        (G['lon'] == 0) |
        (G['second'] != np.round(G['second']))
    )
    return ~bad_gps


def _slow__read_single_v0_9(filename, startMillis, require_gps = True):
    ## this should only be used as a reference
    ## pre-allocate the arrays (more space than is needed)
    M = np.ndarray([15000,12]) # no more than 14400
    G = np.ndarray([15000,11]) # no more than 14400
    D = np.ndarray([750000,2]) # expected number 7.2e5
    d_index = 0
    m_index = 0
    g_index = 0
    millis = startMillis
    ## open the file for reading
    with open(filename, 'r', newline = '', encoding='ascii', errors = 'ignore') as csvfile:
        lines = csv.reader(csvfile, delimiter = ',', quoting = csv.QUOTE_NONE)
        i = 0
        for line in lines:
            ## determine the line type, and skip if it's not necessary data (e.g. debugging info)
            lineType = line[0][0]
            if not (lineType in ['D', 'G', 'M']):
                continue
            ## remove the line type ID and make into a nice array
            if lineType == 'D':
                line[0] = line[0][1:]
            else:
                line = line[1:]
            line = np.array([float(x) for x in line])
            ## unwrap the millis count (always first element of line)
            millis = _unwrap_millis(line[0], millis)
            line[0] = millis
            ## write the line to its matrix
            if lineType == 'D':
                D[d_index,:] = line
                d_index += 1
            elif lineType == 'M':
                M[m_index,:] = line
                m_index += 1
            elif (lineType == 'G') and _check_gps(line):
                G[g_index,:10] = line
                G[g_index,10] = _make_gps_time(line)
                g_index += 1
    if (g_index == 0) and require_gps:
        raise CorruptRawFileNoGPS(filename)
    #pdb.set_trace()
    ## remove unused space in pre-allocated arrays
    D = D[:d_index,:]
    G = pd.DataFrame(G[:g_index,:], columns = ['msPPS', 'msLag', 'year', 'month', 'day', 'hour', \
                                               'minute', 'second', 'lat', 'lon', 't'])
    M = pd.DataFrame(M[:m_index,:], columns = ['millis', 'batt', 'temp', 'A2', 'A3', \
                                               'maxWriteTime', 'minFifoFree', 'maxFifoUsed', \
                                               'maxOverruns', 'gpsOnFlag', 'unusedStack1',\
                                               'unusedStackIdle'])
    ## process data (version-dependent)
    D[:,1] = D[:,1].cumsum()
    return {'data': D, 'metadata': M, 'gps': G}

def _read_several_v0_9(fnList):
    ## initialize the output variables
    G = pd.DataFrame(columns = ['msPPS', 'msLag', 'year', 'month', 'day', 'hour', 'minute', \
                                'second', 'lat', 'lon', 't'])
    M = pd.DataFrame(columns = ['millis', 'batt', 'temp', 'A2', 'A3', 'maxWriteTime', \
                                'minFifoFree', 'maxFifoUsed', 'maxOverruns', 'gpsOnFlag', \
                                'unusedStack1', 'unusedStackIdle'])
    D = np.ndarray([0,2]) # expected number 7.2e5
    num_filler = np.zeros(len(fnList))
    header = pd.DataFrame.from_dict({'file': fnList,
                                     'SN':['' for fn in fnList],
                                     'lat': num_filler,
                                     'lon': num_filler,
                                     'start_ms': num_filler,
                                     'end_ms': num_filler,
                                     'drift_slope': num_filler,
                                     'drift_intercept': num_filler,
                                     'drift_slope_stderr': num_filler,
                                     'drift_resid_std': num_filler,
                                     'drift_resid_MAD': num_filler,
                                     'num_gps_pts': num_filler
                                     })
    ## loop through the files
    startMillis = 0
    for i,fn in enumerate(fnList):
        print('File ' + str(i+1) + ' of ' + str(len(fnList)) + ': ' + fn)
        try:
            L = _read_single_v0_9(fn, startMillis)
        except CorruptRawFileNoGPS:
            print('No GPS data in ' + fn + ', skipping')
        except:
            print('Failed to read ' + fn + ', skipping')
            _breakpoint()
        else:
            if(L['data'][0,0] < startMillis):
                L['metadata'].millis += 2**13
                L['gps'].msPPS += 2**13
                L['data'][:,0] += 2**13
            M = pd.concat((M, L['metadata']))
            G = pd.concat((G, L['gps']))
            D = np.vstack((D, L['data']))
            linreg, num_gps_nonoutliers, MAD_nonoutliers = _robust_regress(L['gps'].msPPS, L['gps'].t)
            resid = L['gps'].t - (linreg.intercept + linreg.slope * L['gps'].msPPS)
            startMillis = D[-1,0]
            header.loc[i, 'lat'] = np.median(L['gps']['lat'])
            header.loc[i, 'lon'] = np.median(L['gps']['lon'])
            header.loc[i, 'start_ms'] = L['data'][0,0] # save this as a millis first, then convert
            header.loc[i, 'end_ms'] = L['data'][-1,0]
            header.loc[i, 'SN'] = _read_SN(fn)
            header.loc[i, 'drift_slope'] = linreg.slope
            header.loc[i, 'drift_intercept'] = linreg.intercept
            header.loc[i, 'drift_slope_stderr'] = linreg.stderr
            header.loc[i, 'drift_resid_std'] = np.std(resid)
            header.loc[i, 'drift_resid_MAD'] = np.max(np.abs(resid))
            header.loc[i, 'num_gps_pts'] = len(L['gps'].msPPS)
            header.loc[i, 'drift_resid_MAD_nonoutliers'] = MAD_nonoutliers
            header.loc[i, 'num_gps_nonoutliers'] = num_gps_nonoutliers
        ## end of fn loop
    return {'metadata':M, 'gps':G, 'data': D, 'header': header}

def _robust_regress(x, y, z=2):
    # goal: a linear regression that is robust to RARE outliers, especially for GPS data
    # scipy.stats.theilslopes (median-based) looks problematic because the median is only affected
    # by the central data point and doesn't benefit from the other samples' information. Also, GPS
    # data slopes are weirdly distributed.
    # In this function, z is the z-score (number of standard deviations) for defining outliers.

    ## Calculate regression line and residuals.
    linreg = scipy.stats.linregress(x, y)
    resid = y - (linreg.intercept + x * linreg.slope)

    ## If any are found to be outliers, remove them and recalculate recursively.
    outliers = np.abs(resid) > (z*np.std(resid))
    if any(outliers):
        return _robust_regress(x[~outliers], y[~outliers], z)
    else:
        return (linreg, len(x), np.max(np.abs(resid)))
    
def _assign_times(L):
    fnList = L['header'].file
    if L['gps'].shape[0] == 0:
        raise Exception('No GPS data in files ' + fnList[0] + '-' + fnList[-1] + '; stopping conversion')
    
    G = _reformat_GPS(L['gps'])
    #breakpoint()
    try:
        breaks = _find_breaks_(L)
    except:
        raise CorruptRawFile('Problem between ' + fnList[0] + '-' + fnList[-1] + '; stopping before this interval. Break between recording periods? Corrupt files?')
    #piecewiseTimeFit = _piecewise_regression(np.array(L['gps'].msPPS), np.array(L['gps'].t), breaks)
    piecewiseTimeFit = L['header']
    L['metadata']['t'] = _apply_segments(L['metadata']['millis'], piecewiseTimeFit)
    header = L['header']
    #header.SN = SN
    header['t1'] = _apply_segments(header.start_ms, piecewiseTimeFit)
    header['t2'] = _apply_segments(header.end_ms, piecewiseTimeFit)
    L['header'] = header
    
    ## interpolate data to equal spacing to make obspy trace
    D = L['data']
    D = np.hstack((D, _apply_segments(D[:,0], piecewiseTimeFit).reshape([D.shape[0],1])))
    timing_info = [L['gps'], L['data'], breaks, piecewiseTimeFit]
    L['data'] = _interp_time(D) # returns stream, populates known fields: channel, delta, and starttime
    L['gps'] = G
    return (L, timing_info)

    
def read_gem(nums = np.arange(10000), path = './', SN = '', units = 'Pa', bitweight = np.NaN, bitweight_V = np.NaN, bitweight_Pa = np.NaN, verbose = True, network = '', station = '', location = '', return_debug_output = False):
    """
    Read raw Gem files.

    Parameters
    ----------
    nums : list or np.array of integers
        Numbers of raw Gem files to read. By default, it reads all files
        in 'path' for the specified serial number.
    
    path : str, default '.' 
        Path of folder containing raw Gem data files to read.

    SN : str
        One Gem serial number to read. Use a loop to read multiple Gems.

    units : str, default 'Pa'
        Desired output units. Options are 'Pa' (Pascals), 'V' (Volts), 
        or 'counts'.

    bitweight : float
        The value of each count when converting between counts and other 
        units (typically Pascals, possibly Volts). By default, it looks
        up the correct bitweight given the Gem's serial number and gain 
        configuration. Leave this blank unless the Gem has been modified
        in a way that changes the bitweight.

    bitweight_V : float
        The value of each count when converting between counts and Volts. 
        By default, it looks up the correct bitweight given the Gem's 
        serial number and gain configuration. Leave blank unless the Gem
        has been modified in a way that changes the voltage bitweight.

    bitweight_Pa : float
        The value of each count when converting between counts and 
        Pascals. By default, it looks up the correct bitweight given the 
        Gem's serial number and gain configuration. Leave this blank 
        unless the Gem has been modified in a way that changes the 
        pressure bitweight.

    verbose : boolean, default True
        Whether to print verbose progress messages to the screen.

    network : str
        Two-character name of the sensor network. Leaving this blank is
        usually fine in subsequent data processing.

    station : str
        Name of the station (up to five characters) to assign to the 
        data. If not provided, uses the Gem's serial number as the 
        station ID.

    location : str
        Two-character location code for this Gem. Leaving this blank is 
        usually fine in subsequent data processing.

    Returns
    -------
    dict with keys:

        - data : obspy.Stream, infrasound data
        - header : pandas.dataframe, information on raw Gem data files
        - metadata : pandas.dataframe, state-of-health and other metadata time series
        - gps : pandas.dataframe, GPS timing and location values

    Note
    ----
    All sample times involving the Gem (and most other passive 
    seismic/acoustic data) are in UTC; time zones are not supported.
    """
    
    if(len(station) == 0):
        station = SN
    fnList = _find_nonmissing_files(path, SN, nums)
    if len(fnList) == 0:
        raise MissingRawFiles
    try:
        version = _read_format_version(fnList[0])
        config = _read_config(fnList[0])
    except:
        raise CorruptRawFile(fnList[0])
    if version == '0.9':
        L = _read_several_v0_9(fnList)
    elif version == '0.85C':
        L = _read_several_v0_9(fnList) # same function works for both
    elif (version == '0.85') | (version == '0.8') :
        raise Exception(fnList[0] + ': Obsolete data format ' + version + ' not yet supported')
    else:
        raise Exception(fnList[0] + ': Invalid or missing data format')

    ## stop early if we don't have data to process
    if len(L['data']) == 0:
        return L
    
    L, timing_info = _assign_times(L)
    
    for tr in L['data']:
        ## populate the rest of the trace stats
        tr.stats.station = station
        tr.stats.location = location # this may well be ''
        tr.stats.network = network # can be '' for now and set later
    ## add bitweight and config info to header
    bitweight_info = get_bitweight_info(SN, config)
    header = L['header']
    for key in bitweight_info.keys():
        L['header'][key] = bitweight_info[key]
    for key in config.keys():
        L['header'][key] = config[key]
    L['header']['file_format_version'] = version

    ## done processing
    if return_debug_output:
        L['debug_output'] = timing_info
    return L

ReadGem = read_gem ## alias, v1.0.0




#########################################################
def _interp_time(data, t1 = -np.Inf, t2 = np.Inf):
    eps = 0.001 # this might need some adjusting to prevent short data gaps
    ## break up the data into continuous chunks, then round off the starts to the appropriate unit
    ## t1 is the first output sample; should be first integer second after or including the first sample (ceiling)
    w_nonnan = ~np.isnan(data[:,2])
    t_in = data[w_nonnan,2]
    p_in = data[w_nonnan,1]
    #_breakpoint()
    t1 = np.trunc(t_in[t_in >= t1][0]+1-0.01) ## 2019-09-11
    t2 = t_in[t_in <= (t2 + .01 + eps)][-1] # add a sample because t2 is 1 sample before the hour
    ## R code here had code to catch t2 <= t1. should add that.
    breaks_raw = np.where(np.diff(t_in) > 0.015)[0]
    breaks = breaks_raw[(t_in[breaks_raw] > t1) & (t_in[breaks_raw+1] < t2)]
    starts = np.hstack([t1, t_in[breaks+1]]) # start times of continuous chunks
    ends = np.hstack([t_in[breaks], t2]) # end times of continuous chunks
    w_same = (starts != ends)
    starts = starts[w_same]
    ends = ends[w_same]
    starts_round = np.trunc(starts)
    ends_round = np.trunc(ends+eps+1)
    ## make an output time vector excluding data gaps, rounded to the nearest samples
    #t_interp = np.zeros(0)
    output = obspy.Stream()
    for i in range(len(starts_round)):
        w = (t_in >= starts[i]) & (t_in <= ends[i])
        try:
            f = scipy.interpolate.CubicSpline(t_in[w], p_in[w])
        except:
            _breakpoint()
            if not _debug:
                raise(Exception('_interp_time failed between ' +str(obspy.UTCDateTime(starts[i])) +\
                                ' and ' + str(obspy.UTCDateTime(ends[i]))))
        t_interp = np.arange(starts[i], ends[i] + eps, 0.01)
        p_interp = np.array(f(t_interp).round(), dtype = 'int32')
        tr = obspy.Trace(p_interp)
        tr.stats.starttime = t_interp[0]
        tr.stats.delta = 0.01
        tr.stats.channel = 'HDF'
        output += tr
    return output
## old code in and below the for loop. no good because a trace is only good for a continuous block. We may have a break, so we need a stream.
#        t_interp = np.concatenate([t_interp, t_new])
#    t_interp = t_interp[(t_interp >= (t1-eps)) & (t_interp < (t2 + eps))]
#    ## interpolate to find pressure at these sample times
#    f = scipy.interpolate.CubicSpline(t_in, p_in)
#    p_interp = np.array(f(t_interp).round(), dtype = 'int32')
#    tr = obspy.Trace(p_interp)
#    tr.stats.starttime = t_interp[0]
#    tr.stats.delta = 0.01
#    tr.stats.channel = 'HDF'
#    return tr

###############################################################
#  version bitweight_Pa  bitweight_V min_SN max_SN
#1    0.50  0.003256538 7.362894e-08      3      7
#2    0.70  0.003256538 7.362894e-08      8     14
#3    0.80  0.003256538 7.362894e-08     15     19
#4    0.82  0.003256538 7.362894e-08     20     37
#5    0.90  0.003543324 7.362894e-08     38     40
#6    0.91  0.003543324 7.362894e-08     41     43
#7    0.92  0.003543324 7.362894e-08     44     46
#8    0.98  0.003501200 7.275362e-08     47     49
#9    0.99  0.003501200 7.275362e-08     50    54
#10    0.991  0.003501200 7.275362e-08     52    54
#11    0.992  0.003501200 7.275362e-08     55    57
#12    1.00  0.003501200 7.275362e-08     58    107
#13    1.00  0.003501200 7.275362e-08     108    Inf
#'bitweight_Pa': [0.003256538, 0.003256538, 0.003256538, 0.003256538, 0.003543324, 0.003543324, 0.003543324, 0.003501200, 0.003501200, 0.003501200, 0.003501200, 0.003501200, 0.003501200],

#bitweight_V = [7.362894e-08, 7.362894e-08, 7.362894e-08, 7.362894e-08, 7.362894e-08, 7.362894e-08, 7.362894e-08, 7.275362e-08, 7.275362e-08, 7.275362e-08, 7.275362e-08,7.275362e-08,7.275362e-08]
__AVCC__ = np.array([3.373, 3.373,  3.373,  3.373,  3.373,  3.373,  3.373, 7.275362e-08, 7.275362e-08, 7.275362e-08, 7.275362e-08, 7.275362e-08, 7.275362e-08]),
def __AVCC__(version):
    if version >= 0.90:
        return 3.1
    else:
        return 3.373

def __gain__(version):
    if version >= 0.98:
        return 1 + 50/0.470
    else:
        return 1 + 49.4/0.470

def get_gem_specs(SN):
    """
    Retrieve specs for a given Gem serial number. 
    
    Parameters
    ----------
    SN : str or int
    
    Returns
    -------
    dict with the following elements:

        - version : Gem version number
        - bitweight_V : voltage resolution for this Gem version [Volts per count]
        - bitweight_Pa : pressure resolution for this Gem version [Pascals per count]
    """
    versionTable = {'version': np.array([0.5, 0.7, 0.8, 0.82, 0.9, 0.91, 0.92, 0.98, 0.99, 0.991, 0.992, 1, 1.01]),
                    'min_SN': np.array([3, 8, 15, 20, 38, 41, 44, 47, 50, 52, 55, 58, 108]),
                    'max_SN': np.array([7, 14, 19, 37, 40, 43, 46, 49, 51, 54, 57, 107, np.Inf])
    }
    version = versionTable['version'][(int(SN) >= versionTable['min_SN']) & (int(SN) <= versionTable['max_SN'])][0]
    bitweight_V = 0.256/2**15/__gain__(version)
    sensitivity = __AVCC__(version)/7.0 * 45.13e-6 # 45.13 uV/Pa is with 7V reference from Marcillo et al., 2012
    return { 'version': version,
             'bitweight_V': bitweight_V,
             'bitweight_Pa': bitweight_V/sensitivity
    }


## Eventually, get_bitweight_info should perform all the functions of the R gemlog equivalent...but not yet.
def get_bitweight_info(SN, config, units = 'Pa'):
    if config['adc_range'] == 0: # high gain
        multiplier = 1
    elif config['adc_range'] == 1: # low gain
        multiplier = 2
    else:
        #pdb.set_trace()
        raise BaseException('Invalid Configuration')
    specs = get_gem_specs(SN)
    specs['bitweight_Pa'] *= multiplier
    specs['bitweight_V'] *= multiplier
    if units.upper() == 'PA':
        specs['bitweight'] = specs['bitweight_Pa']
    elif units.upper() == 'V':
        specs['bitweight'] = specs['bitweight_V']
    elif units.lower() == 'counts':
        specs['bitweight'] = 1
    else:
        raise BaseException('Invalid Units')
    return specs

def _reformat_GPS(G_in):
    t = [obspy.UTCDateTime(tt) for tt in G_in['t']]
    date = [tt.julday + tt.hour/24.0 + tt.minute/1440.0 + tt.second/86400.0 for tt in t]
    G_dict = {'year': np.array([int(year) for year in G_in.year]),
              'date': np.array(date),
              'lat': np.array(G_in.lat),
              'lon': np.array(G_in.lon),
              't': np.array(t)}
    return pd.DataFrame.from_dict(G_dict)

def _find_breaks_(L):
    ## breaks are specified as their millis for comparison between GPS and data
    ## sanity check: exclude suspect GPS tags
    t = np.array([obspy.UTCDateTime(tt) for tt in L['gps'].t])
    tPad = np.array([t[0]] + list(t) + [t[-1]])
    try:
        badTags = ((t > obspy.UTCDateTime.now()) | # no future dates 
                   (L['gps'].year < 2015) | # no years before the Gem existed
                   ((np.abs(t - tPad[:-2]) > 86400) & (np.abs(t - tPad[2:]) > 86400)) | # exclude outliers
                   (L['gps'].lat == 0) | # exclude points within ~1m of the equator
                   (L['gps'].lon == 0)) # exclude points within ~1m of the prime meridian
        L['gps'] = L['gps'].iloc[np.where(~badTags)[0],:]
    except:
        _breakpoint()
    tD = np.array(L['data'][:,0])
    dtD = np.diff(tD)
    starts = np.array([])
    ends = np.array([])
    dataBreaks = np.where((dtD > 25) | (dtD < 0))[0]
    if 0 in dataBreaks:
        tD = tD[1:]
        dtD = dtD[1:]
        dataBreaks = dataBreaks[dataBreaks != 0] - 1
    if (len(dtD) - 1) in dataBreaks:
        tD = tD[:-1]
        dtD = dtD[:-1]
        dataBreaks = dataBreaks[dataBreaks != len(dtD)]
    #_breakpoint()
    for i in dataBreaks:
        starts = np.append(starts, np.max(tD[(i-1):(i+2)]))
        ends = np.append(ends, np.min(tD[(i-1):(i+2)]))
    tG = np.array(L['gps'].t).astype('float')
    mG = np.array(L['gps'].msPPS).astype('float')
    dmG_dtG = np.diff(mG)/np.diff(tG)
    gpsBreaks = np.argwhere(np.isnan(dmG_dtG) | # missing data...unlikely
                         ((np.diff(tG) > 50) & ((dmG_dtG > 1000.1) | (dmG_dtG < 999.9))) | # drift between cycles
                         ((np.diff(tG) <= 50) & ((dmG_dtG > 1002) | (dmG_dtG < 998))) # most likely: jumps within a cycle
                        )
    for i in gpsBreaks:
        i = int(i)
        ## This part is tricky: what if a gpsEnd happens between a dataEnd and dataStart?
        ## Let's be conservative: if either the gpsEnd or gpsStart is within a bad data interval, or what if they bracket a bad data interval?
        ## choose them so that they exclude the most data
        #breakpoint()
        overlaps = (mG[i] <= starts) & (mG[i+1] >= ends)
        if(np.any(overlaps)):
            w = np.argwhere(overlaps)
            starts[w] = max([starts[w], mG[(i-1):(i+2)].max()])
            ends[w] = max([ends[w], mG[(i-1):(i+2)].min()])
        else:
            wmin = np.argwhere(tG > tG[i])
            try:
                starts = np.append(starts, mG[wmin][tG[wmin] == tG[wmin].min()])
            except:
                _breakpoint()
            wmax = np.argwhere(tG < tG[i+1])
            ends = np.append(ends, mG[wmax][tG[wmax] == tG[wmax].max()])
    starts = np.append(tD.min(), starts)
    ends = np.append(ends, tD.max())
    return {'starts':starts, 'ends':ends}

def _apply_segments(x, model):
    y = np.zeros(len(x))
    y[:] = np.nan
    for i in range(len(model['start_ms'])):
        w = (x >= model['start_ms'][i]) & (x <= model['end_ms'][i])
        y[w] = model['drift_intercept'][i] + model['drift_slope'][i] * x[w]
    return y

def _piecewise_regression(x, y, breaks):
    output = {'slope': [], 'intercept':[], 'r':[], 'p':[], 'stderr':[], 'start_ms': [], 'end_ms':[]}
    for i in range(len(breaks['starts'])):
        w = np.where((x >= breaks['starts'][i]) & (x <= breaks['ends'][i]))[0]
        if len(w) == 0: # skip this interval if it doesn't contain data
            continue
        l = scipy.stats.linregress(x[w], y[w])
        output['slope'].append(l.slope)
        output['intercept'].append(l.intercept)
        output['r'].append(l.rvalue)
        output['p'].append(l.pvalue)
        output['stderr'].append(l.stderr)
        output['start_ms'].append(breaks['starts'][i])
        output['end_ms'].append(breaks['ends'][i])
    return output
        
#Concatenate Gem files with no GPS data
