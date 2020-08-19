## gemlog changes:
## new gps column (time)
## metadata time format change

## known potential problems:
## missing/repeated samples at block boundaries due to InterpGem at the beginning
## bitweight missing in output mseed
## not starting files on the hour
## doesn't handle nearly-empty raw files well

## fixed issues:
## from gemlog ReadGemv0.85C (and others too): NaNs in L$gps due to unnecessary and harmful doubling of wna indexing. Also, added python code to drop NaNs.
import pdb
import warnings
import numpy as np
from numpy import NaN, Inf
import os, glob, csv, time, scipy
import pandas as pd
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
import obspy
import datetime
import sys

_debug = True

def _breakpoint():
    if _debug: # skip if we aren't in debug mode
        pdb.set_trace()
#####################
def Convert(rawpath = '.', convertedpath = 'converted', metadatapath = 'metadata', metadatafile = '', gpspath = 'gps', gpsfile = '', t1 = -Inf, t2 = Inf, nums = NaN, SN = '', bitweight = NaN, units = 'Pa', time_adjustment = 0, blockdays = 1, fileLength = 3600, station = '', network = '', location = '', fmt = 'MSEED'):
    ## bitweight: leave blank to use default (considering Gem version, config, and units). This is preferred when using a standard Gem (R_g = 470 ohms)
    
    ## make sure the raw directory exists and has real data
    assert os.path.isdir(rawpath), 'Raw directory ' + rawpath + ' does not exist'
    assert len(glob.glob(rawpath + '/FILE' +'[0-9]'*4 + '.???')) > 0, 'No data files found in directory ' + rawpath

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
    nums.sort()

    ## start at the first file in 'nums'
    n1 = np.min(nums)
  
    ## read the first set of up to (24*blockdays) files
    L = NewGemVar()
    while((L['data'].count() == 0) & (n1 <= max(nums))): ## read sets of files until we get one that isn't empty
        nums_block = nums[(nums >= n1) & (nums < (n1 + (12*blockdays)))] # files are 2 hours, so 12 files is 24 hours
        L = ReadGem(nums_block, rawpath, SN = SN, network = network, station = station, location = location)
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
        gpsfile = makefilename(gpspath, SN, 'gps')

  
    if(len(metadatafile) == 0):
        if(not os.path.isdir(metadatapath)):
            try:
                os.makedirs(metadatapath)
            except:
                print('Failed to make directory ' + metadatapath)
                sys.exit(2)
        metadatafile = makefilename(metadatapath, SN, 'metadata')
  
    ## if the converted directory does not exist, make it
    if(not os.path.isdir(convertedpath)):
        try:
            os.makedirs(convertedpath)
        except:
            print('Failed to make directory ' + metadatapath)
            sys.exit(2)
  
    ## start metadata and gps files
    metadata = L['metadata']   
    gps = L['gps']
    metadata.to_csv(metadatafile, index=False) ## change to metadata format. need to make ScanMnetadata compatible with both

    wgps = (gps['t'] > (t1 - 1)) ## see ReadGemPy to get gps.t working. update gemlog accordingly.
    if(len(wgps) > 0):
        gps[wgps].to_csv(gpsfile, index=False)

    writeHour = max(t1, p[0].stats.starttime)
    writeHour = WriteHourMS(p, writeHour, fileLength, bitweight, convertedpath, fmt=fmt)
    
    ## read sets of (12*blockdays) files until all the files are converted
    while(True):
        ## check to see if we're done
        if(n1 > np.max(nums)):# & len(p) == 0):
            break # out of raw data to convert
        if((t1 > t2) & (not np.isnan(t1 > t2))):
            break # already converted the requested data
        ## load new data if necessary
        tt2 = min(t2, truncUTC(t1, 86400*blockdays) + 86400*blockdays)
        #print([tt2, n1, nums, SN])
        while((p[-1].stats.endtime < tt2) & (n1 <= max(nums))):
            L = ReadGem(nums[(nums >= n1) & (nums < (n1 + (12*blockdays)))], rawpath, SN = SN)
            #pdb.set_trace()
            if(len(L['data']) > 0):
                p = p + L['data']
                p.merge()
            #print(p)
            n1 = n1 + (12*blockdays) # increment file counter
            if(len(L['data']) == 0):
                next # skip ahead if there aren't any readable data files here
                
            ## process newly-read data
            if(any(L['header'].SN != SN) | any(L['header'].SN.apply(len) == 0)):
                #_breakpoint()
                w = (L['header'].SN != SN) | (L['header'].SN.apply(len) == 0)
                print('Wrong or missing serial number(s): ' + L['header'].SN[w] + ' : numbers ' + str(nums[np.logical_and(nums >= n1, nums < (n1 + (12*blockdays)))][w]))
            
            ## start metadata and gps files
            metadata = L['metadata']   
            gps = L['gps']
            metadata.to_csv(metadatafile, index=False, mode='a', header=False)
            wgps = (gps['t'] > (t1 - 1)) ## see ReadGemPy to get gps.t working. update gemlog accordingly.
            ## update the gps file
            if(len(wgps) > 0):
                gps.to_csv(gpsfile, index=False, mode='a', header=False)
                
        ## run the conversion and write new converted files
        #if((pp.stats.endtime >= t1) & (pp.stats.starttime <= tt2))):
        while((writeHour + fileLength) <= p[-1].stats.endtime):
            writeHour = WriteHourMS(p, writeHour, fileLength, bitweight, convertedpath, fmt=fmt)
            
        ## update start time to convert
        p.trim(writeHour, t2)
        t1 = truncUTC(tt2+(86400*blockdays) + 1, 86400*blockdays)
    ## while True
    ## done reading new files. write what's left and end.
    while((writeHour <= p[-1].stats.endtime) & (len(p) > 0)):
        writeHour = WriteHourMS(p, writeHour, fileLength, bitweight, convertedpath, fmt=fmt)
        p.trim(writeHour, t2)
        p = p.split()
        if(len(p) > 0):
            writeHour = p[0].stats.starttime
        else:
            break

def WriteHourMS(p, writeHour, fileLength, bitweight, convertedpath, writeHourEnd = np.nan, fmt='mseed'):
    #pdb.set_trace()
    if(np.isnan(writeHourEnd)):
        writeHourEnd = truncUTC(writeHour, fileLength) + fileLength
    pp = p.copy()
    pp.trim(writeHour, writeHourEnd)
    pp = pp.split() ## in case of data gaps ("masked arrays", which fail to write)
    #_breakpoint()
    for tr in pp:
        tr.stats.calib = bitweight
        fn = MakeFilenameMS(tr, fmt)
        if(len(tr) > 0):
            print(tr)
            if(fmt.lower() == 'wav'):
                ## this is supposed to work for uint8, int16, and int32, but actually only works for uint8. obspy bug?
                    tr.data = np.array(tr.data, dtype = 'uint8')# - np.min(tr[i].data)
                    tr.write(convertedpath +'/'+ fn, format = 'WAV', framerate=7000, width=1) 
            else:
                tr.write(convertedpath +'/'+ fn, format = fmt, encoding=10) # encoding 10 is Steim 1
        #mseed_core._write_mseed(pp, convertedpath +'/'+ fn, format = 'MSEED', encoding=10)

    writeHour = writeHourEnd
    return writeHour
    
## DONE
####################################
## test command
#rawpath = '/home/jake/Work/Gem_Tests/2019-05-29_RoofTestIsolation/raw/'
#SN = '051'
#Convert(rawpath = rawpath, SN = SN, nums = range(14, 15)) 
#Convert(rawpath = rawpath, SN = SN, nums = range(6,8))
#4,15; 5,15; 6,8; 8,10; :ValueError: cannot convert float NaN to integer
#10,12: no error

####################################

def truncUTC(x, n=86400):
    return obspy.core.UTCDateTime(int(float(x)/n)*n)#, origin='1970-01-01')

def makefilename(dir, SN, dirtype):
    n = 0
    fn = dir + '/' + SN + dirtype + '_'+ f'{n:03}' + '.txt'
    while(os.path.exists(fn)):
        n = n + 1
        fn = dir + '/' + SN + dirtype + '_' + f'{n:03}' + '.txt'

    return fn


def MakeFilenameMS(pp, fmt):
    t0 = pp.stats.starttime
    return f'{t0.year:04}' + '-' +f'{t0.month:02}' + '-' +f'{t0.day:02}' + 'T' + f'{t0.hour:02}' + ':' + f'{t0.minute:02}' + ':' + f'{t0.second:02}' + '.' + pp.id + '.' + fmt.lower()
#import pdb

def NewGemVar():
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


def MakeDB(path, pattern = '*', savefile = './DB.csv'):
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

def CalcStationStats(DB, t1, t2):
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


## 55 (3.03), 84 (4.37), 108 (2.04), 49 (1.78), others (1.3-1.6)

#L55=gemlog.ReadGemPy(nums=np.arange(6145,6151),SN='055', path = 'raw')

################################################
def ReadSN(fn):
    SN_line = pd.read_csv(fn, delimiter = ',', skiprows = 4, nrows=1, dtype = 'str', names=['s', 'SN'])
    SN = SN_line['SN'][0]
    return SN

def ReadVersion(fn):
    versionLine = pd.read_csv(fn, delimiter = ',', nrows=1, dtype = 'str', names=['s'])
    version = versionLine['s'][0][7:]
    return version
    
def ReadConfig(fn):
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

def fn2nums(fn_list):
    nums = []
    for i, fn in enumerate(fn_list):
        nums[i] = int(fn[-8:-5])
    return nums


def FindRightFiles(path, SN, nums):
    ## list all Gem files in the path
    fnList = glob.glob(path + '/' + 'FILE[0-9][0-9][0-9][0-9].[0-9][0-9][0-9]')
    fnList.sort()
    fnList = np.array(fnList)
    
    ## find out what all the files' SNs are
    ext = np.array([x[-3:] for x in fnList])
    for i in range(len(ext)):
        if ext[i] == 'TXT':
            ext[i] = ReadSN(fnList[i])
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
        if os.path.getsize(fn) > 0:
            goodNonemptyFnList.append(fn)
    if(len(goodNonemptyFnList) == 0):
        ## warning
        print('No non-empty files')
        return []
    if(len(goodNonemptyFnList) < len(goodFnList)):
        print('Some files are empty, skipping')
    return goodNonemptyFnList    



def UnwrapMillis(new, old, maxNegative = 2**12, rollover = 2**13):
    ## maxNegative is the greatest allowable negative difference.
    ## negative differences can happen between data and GPS lines, or between metadata and data lines.
    return old + ((new - (old % rollover) + maxNegative) % rollover) - maxNegative

def CheckGPS(line): # return True if GPS line is good
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


def MakeGPSTime(line):
    try:
        return obspy.UTCDateTime(int(line[2]), int(line[3]), int(line[4]), int(line[5]), int(line[6]), int(line[7]))
    except:
        return np.NaN

def MillisToTime(L):
    G = L['gps']
    D = L['data']
    coefficients = np.polyfit(G.msPPS, G.t, 3)    
    #print(coefficients)
    pf = np.poly1d(coefficients)
    return pf

def ReadGem_v0_9_single(filename, offset=0):
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
    # skiprows is important so that the header doesn't force dtype=='object'
    # the C engine for pd.read_csv is fast but crashes sometimes. Use the python engine as a backup.
    try:
        df = pd.read_csv(filename, names=range(13), low_memory=False, skiprows=6)
    except:
        df = pd.read_csv(filename, names=range(13), engine='python', skiprows=6, error_bad_lines = False, warn_bad_lines = False)
    df['linetype'] = [value[0] for value in df[0]]
    df = df.loc[df.loc[:,'linetype'].isin(['D', 'M', 'G']), :]
    #df = df.loc[df['linetype'].isin(['D', 'M', 'G']), :]
    ## most of the runtime is before here
    # unroll the ms rollover sawtooth
    df['millis-sawtooth'] = np.where(df['linetype'] == 'D',df[0].str[1:],df[1]).astype(int)
    df['millis-stairstep'] = (df['millis-sawtooth'].diff() < -2**12).cumsum()
    df['millis-stairstep'] -= (df['millis-sawtooth'].diff() > 2**12).cumsum()
    df['millis-stairstep'] *= 2**13
    df['millis-corrected'] = df['millis-stairstep'] + df['millis-sawtooth']
    first_millis = df['millis-corrected'].iloc[0]
    df['millis-corrected'] += (offset-first_millis) + ((first_millis-(offset % 2**13)+2**12) % 2**13) - 2**12
    #df['millis-corrected'] += (offset-first_millis) - ((first_millis - offset) % 2**13)

    # groupby is somewhat faster than repeated subsetting like
    # df.loc[df['linetype'] == 'D', :], ...
    grouper = df.groupby('linetype')
    D = grouper.get_group('D')
    G = grouper.get_group('G')
    M = grouper.get_group('M')

    # pick out columns of interest and rename
    D_cols = ['msSamp', 'ADC']
    G_cols = ['msPPS', 'msLag', 'year', 'month', 'day', 'hour', 'minute',
              'second', 'lat', 'lon']
    M_cols = ['millis', 'batt', 'temp', 'A2', 'A3',
              'maxWriteTime', 'minFifoFree', 'maxFifoUsed',
              'maxOverruns', 'gpsOnFlag', 'unusedStack1', 'unusedStackIdle']

    # column names are currently integers (except for the calculated cols)
    D = D[['millis-corrected', 1]]
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

    return {'data': np.array(D), 'metadata': M.reset_index().astype('float'), 'gps': G.reset_index().astype('float')}


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


def slow_ReadGem_v0_9_single(fn, startMillis):
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
    with open(fn, 'r', newline = '', encoding='ascii', errors = 'ignore') as csvfile:
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
            millis = UnwrapMillis(line[0], millis)
            line[0] = millis
            ## write the line to its matrix
            if lineType == 'D':
                D[d_index,:] = line
                d_index += 1
            elif lineType == 'M':
                M[m_index,:] = line
                m_index += 1
            elif (lineType == 'G') and CheckGPS(line):
                G[g_index,:10] = line
                G[g_index,10] = MakeGPSTime(line)
                g_index += 1
    #pdb.set_trace()
    ## remove unused space in pre-allocated arrays
    D = D[:d_index,:]
    G = pd.DataFrame(G[:g_index,:], columns = ['msPPS', 'msLag', 'year', 'month', 'day', 'hour', 'minute', 'second', 'lat', 'lon', 't'])
    M = pd.DataFrame(M[:m_index,:], columns = ['millis', 'batt', 'temp', 'A2', 'A3', \
                                               'maxWriteTime', 'minFifoFree', 'maxFifoUsed', \
                                               'maxOverruns', 'gpsOnFlag', 'unusedStack1',\
                                               'unusedStackIdle'])
    ## process data (version-dependent)
    D[:,1] = D[:,1].cumsum()
    return {'data': D, 'metadata': M, 'gps': G}

def ReadGem_v0_9(fnList):
    ## initialize the output variables
    G = pd.DataFrame(columns = ['msPPS', 'msLag', 'year', 'month', 'day', 'hour', 'minute', 'second', 'lat', 'lon', 't'])
    M = pd.DataFrame(columns = ['millis', 'batt', 'temp', 'A2', 'A3', \
                                               'maxWriteTime', 'minFifoFree', 'maxFifoUsed', \
                                               'maxOverruns', 'gpsOnFlag', 'unusedStack1',\
                                               'unusedStackIdle'])
    D = np.ndarray([0,2]) # expected number 7.2e5
    num_filler = np.zeros(len(fnList))
    header = pd.DataFrame.from_dict({'file': fnList,
                                     'SN':['' for fn in fnList],
                                     'lat': num_filler,
                                     'lon': num_filler,
                                     't1': num_filler,
                                     't2': num_filler
                                     })
    ## loop through the files
    startMillis = 0
    for i,fn in enumerate(fnList):
        print('File ' + str(i+1) + ' of ' + str(len(fnList)) + ': ' + fn)
        try:
            L = ReadGem_v0_9_single(fn, startMillis)
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
            startMillis = D[-1,0]
            header.loc[i,'lat'] = np.median(L['gps']['lat'])
            header.loc[i,'lon'] = np.median(L['gps']['lon'])
            header.loc[i, 't1'] = L['data'][0,0] # save this as a millis first, then convert
            header.loc[i,'t2'] = L['data'][-1,0]
        ## end of fn loop
    return {'metadata':M, 'gps':G, 'data': D, 'header': header}


def ReadGem(nums = np.arange(10000), path = './', SN = '', units = 'Pa', bitweight = np.NaN, bitweight_V = np.NaN, bitweight_Pa = np.NaN, verbose = True, network = '', station = '', location = ''):
    if(len(station) == 0):
        station = SN
    ## add asserts, especially for SN
    fnList = FindRightFiles(path, SN, nums)
    version = ReadVersion(fnList[0])
    config = ReadConfig(fnList[0])
    if version == '0.9':
        L = ReadGem_v0_9(fnList)
    elif version == '0.85C':
        L = ReadGem_v0_9(fnList) # will the same function work for both?
    elif (version == '0.85') | (version == '0.8') :
        raise Exception(fnList[0] + ': Obsolete data format ' + version + ' not yet supported')
    else:
        raise Exception(fnList[0] + ': Invalid or missing data format')
    assert L['gps'].shape[0] > 0, 'No GPS data in files ' + fnList[0] + '-' + fnList[-1] + '; stopping conversion'
    M = L['metadata']
    D = L['data']
    G = ReformatGPS(L['gps'])
    #breakpoint()
    breaks = FindBreaks(L)
    piecewiseTimeFit = PiecewiseRegression(np.array(L['gps'].msPPS), np.array(L['gps'].t), breaks)
    M['t'] = ApplySegments(M['millis'], piecewiseTimeFit)
    header = L['header']
    header.SN = SN
    header.t1 = ApplySegments(header.t1, piecewiseTimeFit)
    header.t2 = ApplySegments(header.t2, piecewiseTimeFit)
    D = np.hstack((D, ApplySegments(D[:,0], piecewiseTimeFit).reshape([D.shape[0],1])))
    
    ## interpolate data to equal spacing to make obspy trace
    #_breakpoint()
    st = InterpTime(D) # populates known fields: channel, delta, and starttime
    for tr in st:
        ## populate the rest of the trace stats
        tr.stats.station = station
        tr.stats.location = location # this may well be ''
        tr.stats.network = network # can be '' for now and set later
    ## add bitweight and config info to header
    bitweight_info = GetBitweightInfo(SN, config)
    for key in bitweight_info.keys():
        header[key] = bitweight_info[key]
    for key in config.keys():
        header[key] = config[key]
    header['file_format_version'] = version
    return {'data': st, 'metadata': M, 'gps': G, 'header' : header}






#########################################################
def InterpTime(data, t1 = -np.Inf, t2 = np.Inf):
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

def GemSpecs(SN):
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


## Eventually, GetBitweightInfo should perform all the functions of the R gemlog equivalent...but not yet.
def GetBitweightInfo(SN, config, units = 'Pa'):
    if config['adc_range'] == 0: # high gain
        multiplier = 1
    elif config['adc_range'] == 1: # low gain
        multiplier = 2
    else:
        #pdb.set_trace()
        raise BaseException('Invalid Configuration')
    specs = GemSpecs(SN)
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

def ReformatGPS(G_in):
    t = [obspy.UTCDateTime(tt) for tt in G_in['t']]
    date = [tt.julday + tt.hour/24.0 + tt.minute/1440.0 + tt.second/86400.0 for tt in t]
    G_dict = {'year': np.array([int(year) for year in G_in.year]),
              'date': np.array(date),
              'lat': np.array(G_in.lat),
              'lon': np.array(G_in.lon),
              't': np.array(t)}
    return pd.DataFrame.from_dict(G_dict)


def FindBreaks(L):
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

def ApplySegments(x, model):
    y = np.zeros(len(x))
    y[:] = np.nan
    for i in range(len(model['start'])):
        w = (x >= model['start'][i]) & (x <= model['end'][i])
        y[w] = model['intercept'][i] + model['slope'][i] * x[w]
    return y

def PiecewiseRegression(x, y, breaks):
    output = {'slope': [], 'intercept':[], 'r':[], 'p':[], 'stderr':[], 'start': [], 'end':[]}
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
        output['start'].append(breaks['starts'][i])
        output['end'].append(breaks['ends'][i])
    return output
        
