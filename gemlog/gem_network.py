import numpy as np
import pandas as pd
import glob, obspy, os, warnings, gemlog, sys, argparse, pathlib
#from obspy.clients.nrl import NRL
#from contextlib import contextmanager,redirect_stderr,redirect_stdout
#from os import devnull
#nrl = NRL()

#response = nrl.get_response(sensor_keys = ['Gem', 'Gem Infrasound Sensor v1.0'], datalogger_keys = ['Gem', 'Gem Infrasound Logger v1.0', '0 - 128000 counts/V'])

#@contextmanager
#def _suppress_stdout_stderr():
#    """A context manager that redirects stdout and stderr to devnull"""
#    with open(devnull, 'w') as fnull:
#        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
#            yield (err, out)

def deconvolve_gem_response(data, gain = 'high', sensor_file = '', logger_file = ''):
    """
    Remove the Gem's instrument response
    
    Parameters:
    -----------
    data : either obspy.Stream or obspy.Trace
    Input data to remove response from. Note that if a stream is provided, the Gem response will be 
    deconvolved from all traces, even if some traces were not recorded by Gems!

    gain : str, default 'high'
    If a configuration file was used to set the Gem's programmable gain to half-gain, use 'low'.

    Returns:
    --------
    The input data (either trace or stream) with response removed
    
    Example:
    --------
    import obspy, gemlog
    ## read sample data--this isn't from a Gem, but it works
    tr = obspy.read()[0] 
    tr_deconvolved = gemlog.deconvolve_gem_response(tr)
    """
    ## obspy's tool to read instrument responses always triggers a warning for infrasound sensors
    ## use the warnings package to suppress warnings in this function only
    if type(data) == obspy.Stream:
        for tr in data:
            tr.stats.response = get_gem_response(gain, sensor_file, logger_file)
    elif type(data) == obspy.Trace:
        data.stats.response = get_gem_response(gain, sensor_file, logger_file)
    else:
        raise TypeError(f'data is type {type(data)}; must be obspy.Stream or obspy.Trace')

    ## By now, the response has been attached to all traces. Remove the response and return the data.
    data.remove_response()
    return data

## Create response variables. Do this once on package import so it doesn't have to be run repeatedly.
## obspy's tool to read instrument responses always triggers a warning for infrasound sensors
## use the warnings package to suppress warnings in this function only
#with warnings.catch_warnings():
#    warnings.simplefilter("ignore")
#_response_high_gain_1_0 = nrl.get_response(sensor_keys = ['Gem', 'Gem Infrasound Sensor v1.0'],
#                                           datalogger_keys = ['Gem', 'Gem Infrasound Logger v1.0',
#                                                              '0 - 128000 counts/V']) # may cause warning--ok to ignore
#_response_low_gain_1_0 =  nrl.get_response(sensor_keys = ['Gem', 'Gem Infrasound Sensor v1.0'],
#                                           datalogger_keys = ['Gem', 'Gem Infrasound Logger v1.0',
#                                                              '1 - 64000 counts/V']) # may cause warning--ok to ignore

def _read_response(filename):
    return obspy.read_inventory(filename)[0][0][0].response

def get_gem_response(gain = 'high', sensor_file = '', logger_file = ''):
    """
    Return the Gem's instrument response
    
    Parameters:
    -----------
    gain : str, default 'high'
    If a configuration file was used to set the Gem's programmable gain to half-gain, use 'low'.
    High gain is normal.

    sensor_file : str
    Not normally used; allows user to select a non-standard response (uncommon)

    logger_file : str
    Not normally used; allows user to select a non-standard response (uncommon)

    Returns:
    --------
    The Gem's instrument response as obspy.core.inventory.response.Response

    Example:
    --------
    import obspy, gemlog
    ## read sample data--this isn't from a Gem, but it works as a demo
    tr = obspy.read()[0] 

    ## note that the remainder of this example can be done more directly using 
    gemlog.deconvolve_gem_response

    ## find the Gem's response
    response = gemlog.get_gem_response()

    ## attach the response to the trace. 
    tr.stats.response = response

    ## remove the response from the trace in place
    tr.remove_response()
    
    """
    response_path = os.path.join(os.path.dirname(gemlog.__file__), 'data', 'response')
    if sensor_file == '':
        sensor_file = 'RESP.XX.IS025..BDF.GEMV1.26.0_0022'
    sensor_resp = _read_response(os.path.join(response_path, 'sensor', sensor_file))

    if logger_file == '':
        if gain.lower() == 'high':
            logger_file = 'RESP.XX.GM002..HHZ.GEMINFRAV1.0.100'
        elif gain.lower() == 'low':
            logger_file = 'RESP.XX.GM002..HHZ.GEMINFRAV1.0.100'
        else:
            raise ValueError(f'invalid gain: {gain}')
    response = _read_response(os.path.join(response_path, 'datalogger', logger_file))

    ## code copied from obspy's nrl.get_response() to merge sensor and logger responses
    response.response_stages.pop(0)
    sensor_stage0 = sensor_resp.response_stages[0]
    response.response_stages.insert(0, sensor_stage0)
    response.instrument_sensitivity.input_units = sensor_stage0.input_units
    response.instrument_sensitivity.input_units_description = sensor_stage0.input_units_description
    
    ## Obspy's get_response method fails to set the overall sensitivity correctly for infrasound
    ## responses, so we have to do this manually. If this isn't done, we'll get a warning later.
    ## We have to do this at freq 0.05 Hz for compatibility with the datalogger's nominal gain freq.
    filter_gain_0_05 = np.abs(1-1/(1+0.05j/0.039))
    response.instrument_sensitivity.value = response.response_stages[0].stage_gain * response.response_stages[2].stage_gain * filter_gain_0_05
    return response

## Two issues with obspy warnings with instrument responses:
## issue 1: verifying the response with nrl.get_response. After the sensor and datalogger responses are combined, it runs dl_resp.recalculate_overall_sensitivity() as a sanity check. Unfortunately, this method is unaware of pressure and pascals, and fails if given any option other than displacement, velocity, or acceleration. This can be suppressed with the warnings package, or avoided by just not using nrl.get_response().
## issue 2: removing the response, which ultimately calls response._call_eval_resp_for_frequencies(). This calls compiled C code, and the C code complains because it can't verify the sensitivity. This cannot be suppressed by either the warnings package or the io suppression. However, it can be suppressed by MANUALLY changing the overall instrument gain for the gem response, which is now done at the end of get_gem_response. Warnings should be expected for any other means of getting the gem response.
def _fix_station_info_keys(d):
    ## after reading a station_info file, standardize the key capitalization
    ## this saves the user from having to be picky about capitalizing the header right
    allowed_keys = ['SN', 'network', 'station', 'location', 'elevation']
    for key in allowed_keys:
        ## for the given allowed key, find the existing key that matches it
        w = [test_key.lower()[:2] == key.lower()[:2] for test_key in d.keys()]
        if any(w):
            ## rename that key to the standard name
            d[key] = d.pop(list(d.keys())[np.where(w)[0][0]])
    return d

def _get_station_info(station_info):
    required_keys = ['sn', 'network', 'station', 'location']
    if type(station_info) == str:
        print('Reading file %s' % station_info)
        header_df = pd.read_csv(station_info, nrows = 1, header = None, index_col = False)
        header_list = [header_df[key][0] for key in header_df.keys()]
        if all([i.lower() in (required_keys + ['elevation']) for i in header_list]): # file has header line
            print('File has valid header, using that for column names')
            station_info = pd.read_csv(station_info, dtype = {key:'str' for key in header_list}, keep_default_na = False)
        else:
            if len(header_list) == 4:
                print('File does not have a valid header, using default columns [SN, network, station, location]')
                station_info = pd.read_csv(station_info, names = required_keys, dtype = {key:'str' for key in header_list}, keep_default_na = False, index_col = False)
            elif len(header_list) == 5:
                print('File does not have a valid header, using default columns [SN, network, station, location, elevation]')
                station_info = pd.read_csv(station_info, names = required_keys + ['elevation'], dtype = {key:'str' for key in header_list}, keep_default_na = False, index_col = False)
            else:
                raise Exception('invalid station_info file; must have 4 or 5 columns or valid header')
        ## if any keys are capitalized or abbreviated wrong, correct them
        station_info = _fix_station_info_keys(station_info)
    elif (type(station_info) is not pd.DataFrame):
        raise Exception('invalid station_info input')
    else:
        #station_info = _fix_station_info_keys(station_info)
        if any([key not in station_info.keys() for key in required_keys]):
            raise Exception('invalid station_info input')
    # if location and network fields are blank in file, they are interpreted as NaN and must be
    # turned back into blank
    if 'elevation' not in station_info.keys():
        station_info['elevation'] = np.zeros(station_info.shape[0]) - 9999
        
    if station_info['elevation'].dtype == 'str':
        station_info['elevation'] = station_info['elevation'].astype('float')
        
    station_info.loc[station_info.location.isna(), 'location'] = ''
    station_info.loc[station_info.network.isna(), 'network'] = ''
    station_info.loc[station_info.elevation.isna(), 'elevation'] = -9999
    return station_info

def make_gem_inventory(station_info, coords, response = 'default'):
    """
    Create a station inventory for a Gem dataset.

    Parameters
    ----------
    station_info: str or pandas.DataFrame
        file that contains a table with station definitions (columns for serial number, network,
        station, location, and channel)
    
    coords : pandas.DataFrame
        output of summarize_gps(), with 'elevation' column added

    response: str
        instrument response information (currently only 'default' is supported)

    Returns
    -------
    obspy.Inventory containing station metadata for the dataset
        
    """
    ## ensure that the provided station_info has all the necessary info, and fail if it doesn't
    station_info = _get_station_info(station_info)
    ## ensure that the coords info has all the necessary info, and fail if it doesn't
    if (type(coords) is not pd.DataFrame):
        raise Exception('invalid coords; must be pandas.DataFrame')
    elif ('SN' not in coords.keys()) and any([key not in coords.keys() for key in ['network', 'station', 'location']]):
        raise Exception("invalid coords; must contain key 'SN', or keys 'network', 'station', 'location'")
    elif any([key not in coords.keys() for key in ['lat', 'lon', 'elevation']]):
        raise Exception("invalid coords; must contain keys 'lat', 'lon', 'elevation'")
    
    ## create the inventory and loop through all the networks, stations, and locations in it
    inventory = obspy.Inventory()
    all_networks = _unique(station_info['network'])
    for network_name in all_networks:
        network = obspy.core.inventory.Network(network_name)
        ## loop through all the stations in the current network
        all_stations = _unique(station_info[station_info['network'] == network_name]['station'])
        for station_name in all_stations:
            station = obspy.core.inventory.Station(station_name, latitude = 0, longitude = 0, elevation = 0)
            station.site.name = station_name
            ## loop through all the locations/channels in the current station
            all_locations = _unique(station_info[(station_info['network'] == network_name) & (station_info['station'] == station_name)]['location'])
            for location_name in all_locations:
                ## find the serial number for this location
                SN = station_info[(station_info['network'] == network_name) & (station_info['station'] == station_name) & (station_info['location'] == location_name)]['SN'].iloc[0]
                ## find the coords line corresponding to this location
                if ('SN' in coords.keys()) and any(key not in coords.keys() for key in ['network', 'station', 'location']):
                    line = coords[coords['SN'].astype(int) == int(SN)] # int just in case there's some type discrepancy
                else:
                    #line = coords[(coords['network'] == network_name) & (coords['station'] == station_name) & (coords['location'].astype(int) == int(location_name))] ## this crashes when location is ''
                    line = coords[(coords['network'] == network_name) & (coords['station'] == station_name) & (coords['location'] == location_name)]
                ## if no GPS info for this station is found, skip it
                if line.shape[0] == 0:
                    print('No coords found for %s.%s.%s, skipping' % (network_name, station_name, location_name))
                    continue
                ## We need to extract the coordinate and times for this location.
                lat = line['lat'].iloc[0]
                lon = line['lon'].iloc[0]
                elevation = line['elevation'].iloc[0]
                ## Check to see if the coords include start/end times. If so, pad
                ## them 1 hour to be safe. If not, assume the station is eternal.
                if all([key in coords.keys() for key in ['starttime', 'endtime']]):
                    t1 = obspy.UTCDateTime(line['starttime'].iloc[0]) - 3600 
                    t2 = obspy.UTCDateTime(line['endtime'].iloc[0]) + 3600
                else:
                    t1 = obspy.UTCDateTime('1970-01-01T00:00:00')
                    t2 = obspy.UTCDateTime('9999-12-31T23:59:59')
                ## Assume it's a gem v1.0--safe for now
                equipment = obspy.core.inventory.util.Equipment(serial_number = SN, model = 'Gem Infrasound Logger v1.0', description = 'Gem 1.0 (Infrasound), 0.039-27.1 Hz, 0.0035012 Pa/count')
                channel = obspy.core.inventory.Channel('HDF', location_code = location_name, latitude = lat, longitude = lon, elevation = elevation, depth = 0, response = response, equipments = equipment, start_date = t1, end_date = t2, sample_rate = 100, clock_drift_in_seconds_per_sample = 0, types = ['GEOPHYSICAL'], sensor = equipment, azimuth = 0, dip = 0)
                station.channels.append(channel)
            ## move on to the next station if this one is empty
            if len(station.channels) == 0:
                continue
            ## calculate the current station's coordinate as the mean of its locations,
            ## then append it to the network
            station.latitude = np.mean([channel.latitude for channel in station.channels])
            station.longitude = np.mean([channel.longitude for channel in station.channels])
            station.elevation = np.mean([channel.elevation for channel in station.channels])
            station.start_date = np.min([channel.start_date for channel in station.channels])
            station.end_date = np.max([channel.end_date for channel in station.channels])
            network.stations.append(station)
        ## append the current network to the inventory
        network.start_date = np.min([station.start_date for station in network.stations])
        network.end_date = np.max([station.end_date for station in network.stations])
        inventory.networks.append(network)
    return inventory
    
def rename_files(infile_pattern, station_info, output_dir, output_format = 'mseed', outfile_pattern =
                 '{year}-{mon}-{day}T{hour}_{min}_{sec}.{net}.{sta}.{loc}.{chan}.{fmt}'):
    """
    Rename a set of converted data files, assigning the correct network, station, and location 
    codes (instead of the original code with empty location/network and the station code being the
    serial number).

    Parameters
    ----------
    infile_pattern : str
        glob-type pattern defining the input files

    station_info : str or pandas.DataFrame
        file name for table assigning serial numbers to network, station, and location codes

    output_dir : str
        folder where output files should be written

    output_format : str
        default 'mseed'; 'sac' also works, as do other obspy-supported formats

    outfile_pattern : str
        format of output file names; note the abbreviations, and that 'jd' means 'day of year'

    Returns
    -------
    pandas.DataFrame containing the station_info table.
    """
    station_info = _get_station_info(station_info)
    
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir) # makedirs vs mkdir means if gpspath = dir1/dir2, and dir1 doesn't exist, that dir1 will be created and then di
    infiles = glob.glob(infile_pattern)
    infiles.sort()
    assert len(infiles) > 0, 'No input files provided'
    for i, infile in enumerate(infiles):
        st = obspy.read(infile)
        fileParts = infile.split('/')[-1].split('.')
        SN = fileParts[2]
        try:
            w = np.where(station_info.SN == SN)[0][0]
        except:
            print('skipping ' + infile)
            continue
        network = station_info.network[w]
        station = station_info.station[w]
        location = station_info.location[w]
        t1 = st[0].stats.starttime
        for tr in st:
            tr.stats.network = network
            tr.stats.station = station
            tr.stats.location = location
            t1 = min(t1, tr.stats.starttime)
        #outputFile = output_dir + '/' + '%s.%s.%s.%s.HDF.%s' % (fileParts[0], network, station, location, output_format)
        trace_info_dict = {'year':t1.year, 'mon':t1.month, 'day':t1.day,
                           'jd':t1.julday, 'hour':t1.hour, 'min':t1.minute,
                           'sec':t1.second, 'net':network, 'sta':station,
                           'loc':location, 'chan':'HDF', 'fmt':output_format}
        outfile_pattern = _fix_file_name_digits(outfile_pattern)
        output_file = outfile_pattern.format(**trace_info_dict)
        st.write(output_dir + '/' + output_file, format = output_format)
        print(str(i+1) + ' of ' + str(len(infiles)) + ': ' + infile + ', ' + output_file)
    return station_info

def merge_files_day(infile_path, infile_pattern = '*', outfile_dir = 'merge_file_output'):
    if not os.path.isdir(outfile_dir):
        os.makedirs(outfile_dir) # makedirs vs mkdir means if gpspath = dir1/dir2, and dir1 doesn't exist, that dir1 will be created and then dir2
    infiles = glob.glob(infile_path + '/' + infile_pattern)
    infiles.sort()
    ## find unique year, month, day lists for all files
    cuts = [infile.split('/')[-1].split('.')[0].split('T')[0] + ',' + '.'.join(infile.split('/')[-1].split('.')[1:])
            for infile in infiles]
    cuts = _unique(cuts)
    days = [cut.split(',')[0] for cut in cuts]
    suffixes = [cut.split(',')[1] for cut in cuts]
    
    ## for each item, read all the files and 
    for day, suffix in zip(days, suffixes):
        print(day + ' ' + suffix)
        x = obspy.read(infile_path + '/' + day + '*' + suffix)
        x.merge(fill_value = 'latest', method = 1)
        x.write(outfile_dir + '/' + day + 'T00_00_00.' + suffix)
    return

def _fix_file_name_digits(fn):
    fn = fn.replace('{year}', '{year:04d}')
    fn = fn.replace('{mon}', '{mon:02d}')
    fn = fn.replace('{day}', '{day:02d}')
    fn = fn.replace('{hour}', '{hour:02d}')
    fn = fn.replace('{min}', '{min:02d}')
    fn = fn.replace('{sec}', '{sec:02d}')
    fn = fn.replace('{jd}', '{jd:03d}')
    return fn


# function to get unique values from list while preserving order (set-based shortcut doesn't do this)
def _unique(list1): 
    unique_list = []
    for x in list1:
        # check if exists in unique_list or not 
        if x not in unique_list: 
            unique_list.append(x) 
    return unique_list

## function to exclude outliers by repeatedly calculating standard dev and tossing points outside N standard devs, until none are left
## this matters because occasionally a dataset will start or end with short recordings made elsewhere, which must be excluded from the calculation of station coords
def _remove_outliers(x, N = 5):
    w = np.where((np.abs(x.lat - np.median(x.lat)) < (N * (1e-5+np.std(x.lat)))) & (np.abs(x.lon - np.median(x.lon)) < (N * (1e-5+np.std(x.lon)))))[0]
    if len(w) < x.shape[0]:
        return _remove_outliers(x.iloc[w,:], N=N)
    else:
        return(x)

def read_gps(gps_dir_pattern, SN):
    """
    Read the most up-to-date GPS file for a given serial number.

    Parameters
    ----------
    gps_dir_pattern : str
        Path to the folder containing GPS data, or a glob-style pattern describing multiple folders.

    SN : str
        Serial number of the Gem being examined.

    Returns
    -------
    pandas.DataFrame containing columns year, date (day of year), lat, lon (all floats), and column
    t (obspy.UTCDateTime).
    """
    if type(gps_dir_pattern) is str:
        gpsDirList = sorted(glob.glob(gps_dir_pattern))
    elif type(gps_dir_pattern) is list:
        gpsDirList = gps_dir_pattern
    else:
        raise Exception('invalid gps_dir_pattern type; must be list or str')
    gpsTable = pd.DataFrame(columns=['year', 'date', 'lat', 'lon', 't'])
    for gpsDir in gpsDirList:
        if os.path.isdir(gpsDir):
            fnList = sorted(glob.glob(gpsDir + '/' + SN + '*'))
        else:
            # in case the user provides the gps filename directly. requires that the filename start with SN.
            if pathlib.Path(gpsDir).name[:3] == SN:
                fnList = [gpsDir]
            else:
                fnList = []
        if len(fnList) > 0: # if any gps files matching SN are found, read and append the last
            gpsTable = pd.concat([gpsTable,
                                  pd.read_csv(fnList[-1])], ignore_index=True)
    return gpsTable
ReadLoggerGPS = read_gps # alias; v1.0.0

def summarize_gps(gps_dir_pattern, station_info = None, output_file = None, t1 = None, t2 = None, include_SN = None, exclude_SN = None):
    """
    Read up-to-date GPS data from all Gems in a project, estimate their locations using a robust
    trimmed-mean method.

    Parameters
    ----------
    gps_dir_pattern : str
        Path to folder or glob-style pattern describing folders containing GPS data to review.

    station_info : str
        Path to text file containing table assigning serial numbers to network, station, and 
        location codes.

    output_file : str
        Path to file where output should be written (optional)

    t1 : obspy.UTCDateTime or str
        Start time; ignore gps data before t1 (e.g., '2022-12-31' or '2022-12-31_05:45:00')

    t2 : obspy.UTCDateTime or str
        End time; ignore gps data after t2

    include_SN: list
        Serial numbers to process (default all)

    exclude_SN: list
        Serial numbers not to process (default none)

    Returns
    -------
    pandas.DataFrame containing the following columns:

        - SN: serial number (str)
        - lat: calculated average latitude (float)
        - lon: calculated average longitude (float)
        - lat_SE: standard error of average latitude (float)
        - lon_SE: standard error of average longitude (float)
        - elevation: elevation provided by user (-9999 if not provided) (float)
        - starttime: time of first GPS data (obspy.UTCDateTime)
        - endtime: time of last GPS data (obspy.UTCDateTime)
        - num_samples: number of GPS strings recorded 
    

    If station_info is provided, then the following columns are also included:

        - network: network code (str)
        - station: station code (str)
        - location: location code (str)
        - elevation: elevation provided in station_info (str)
    """
    if t1 is not None:
        t1 = obspy.UTCDateTime(t1)
    else:
        t1 = obspy.UTCDateTime('1970-01-01')
    if t2 is not None:
        t2 = obspy.UTCDateTime(t2)
    else:
        t2 = obspy.UTCDateTime('9999-12-31')

    if type(gps_dir_pattern) is str:
        gpsDirList = sorted(glob.glob(gps_dir_pattern))
        gpsFileList = []
        for gpsDir in gpsDirList:
            gpsFileList += glob.glob(gpsDir + '/' + '*gps*txt')
    elif type(gps_dir_pattern) is list:
        gpsFileList = []
        for fn in gps_dir_pattern:
            if os.path.isdir(fn):
                gpsFileList += glob.glob(fn + '/*gps*txt')
            else:
                gpsFileList.append(fn)
    else:
        raise Exception(f'invalid type for gpsDirList {gpsDirList}; must be str or list')
    
    gpsFileList = sorted(gpsFileList)
    snList = []
    for gpsFile in gpsFileList:
        current_SN = gpsFile.split('/')[-1].split('gps')[0]
        if ((include_SN is None) or (current_SN in include_SN)) and \
           ((exclude_SN is None) or (current_SN not in exclude_SN)):
            snList.append(current_SN)
            
    snList = sorted(_unique(snList))
    coords = pd.DataFrame(columns = ['SN', 'lat', 'lon', 'lat_SE', 'lon_SE', 'starttime', 'endtime', 'num_samples'])
    avgFun = lambda x: np.mean(x)
    seFun = lambda x: np.std(x)/np.sqrt(len(x))
    for i, SN in enumerate(snList):
        print(str(i+1) + ' of ' + str(len(snList)) + ': ' + SN)
        gpsTable = _remove_outliers(read_gps(gps_dir_pattern, SN))
        gpsTable = gpsTable.iloc[np.where((gpsTable.t >= t1) & (gpsTable.t <= t2))[0], :]
        if gpsTable.shape[0] > 0:
            coords = pd.concat([
                coords,
                pd.DataFrame(
                    [[SN, avgFun(gpsTable.lat), avgFun(gpsTable.lon), seFun(gpsTable.lat), seFun(gpsTable.lon), gpsTable.t.min(), gpsTable.t.max(), gpsTable.shape[0]]],
                    columns = ['SN', 'lat', 'lon', 'lat_SE', 'lon_SE', 'starttime', 'endtime', 'num_samples'])
                ], ignore_index = True)
                               
        else:
            print('No non-outliers remaining for gem {SN}')
    if station_info is not None:
        station_info = _get_station_info(station_info)
        network = []
        station = []
        location = []
        elevation = []
        for SN in coords.SN:
            try:
                w = np.where(station_info.SN == SN)[0][0]
                network.append(station_info.network[w])
                station.append(station_info.station[w])
                location.append(station_info.location[w])
                elevation.append(station_info.elevation[w])
            except:
                network.append('')
                station.append('')
                location.append('')
                elevation.append(-9999)
        coords['network'] = network
        coords['station'] = station
        coords['location'] = location
        coords['elevation'] = elevation
    else:
        coords['elevation'] = [-9999 for SN in coords.SN]
            
    if output_file is not None:
        coords.to_csv(output_file, index = False)
    return coords
    
SummarizeAllGPS = summarize_gps # alias, v1.0.0

def summarize_gps_terminal(input = sys.argv[1:]):
    details='''
    The required input gps_folder is the path to the folder where gemconvert saved the gps data.

    t1 and t2 are useful for when the dataset includes data outside the deployment period (e.g., data
    from preliminary or subsequent testing that was left on the disk and processed along with the 
    field data).

    The station_info csv file (optional) defines how the Gem sensors are assigned SEED codes.
    The csv file must including at least the following columns: SN, station. 

    Network and location columns are optional; if not provided, they are assumed to be blank ('').

    Station code must be 1-5 characters; location and network must be 0-2 characters. 
    All Gem recordings are assigned channel code HDF (infrasound sampled at 100 sps).
    Because each Gem is its own data logger, it is recommended to assign each Gem a different station name
    and to leave location codes blank. Additionally, if you plan on submitting your data to IRIS-DMC or
    similar, you should plan to use a network code issued by the DMC for your whole dataset.

    Elevation column (meters) is optional. The available GPS info does not contain elevations, but IRIS-DMC
    requires elevation to be provided; to comply with that requirement, you must determine elevation
    yourself and enter it as a column. If not provided, elevation will be assumed to be -9999.

    For example, if the station_info file consists of the following:
    SN,Network,Station,Location,Elevation
    077,NM,LADR1,,1022
    088,NM,LADR2,,1016

    Then Gem #077 will be assigned SEED code NM.LADR1..HDF, and Gem #088 is NM.LADR2..HDF. You will
    get identical results with the following station_info file that omits location (since they are all blank):
    SN,Network,Station,Elevation
    077,NM,LADR1,1022
    088,NM,LADR2,1016
    
    Finally, the following file is the minimum information that can be processed. Note that this is not
    compliant with IRIS-DMC's requirement for a network code and elevation.
    sn,station
    077,LADR1
    088,LADR2
    '''
    parser = argparse.ArgumentParser(description='Summarize gps information to make ',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog = details)
    parser.add_argument('gps_folder', nargs='+', help='list of gps files (wildcards are allowed, e.g. gps/*)')
    
    parser.add_argument('-o', '--output_file', default='gem_network',
                        help='Name (without file type extension) of output files to write')


    parser.add_argument('-s', '--station_info_file', default=None, help = 'csv file containing additional info on station names (see details below)')

    parser.add_argument('-i', '--include_SNs', help='Serial numbers to include in processing (default: all)')
    parser.add_argument('-x', '--exclude_SNs', help='Serial numbers to exclude from processing (default: none)')

    parser.add_argument('-1', '--t1', default='1970-01-01', help='Time to start processing data (default: beginning of data)')
    parser.add_argument('-2', '--t2', default='9999-12-31', help='Time to stop processing data (default: end of data)')
    args = parser.parse_args(input)
    
    if args.include_SNs is not None and type(args.include_SNs) is str:
        include_SNs = args.include_SNs.split(',')
    else:
        include_SNs = None
    if (args.exclude_SNs is not None) and type(args.exclude_SNs) is str:
        exclude_SNs = args.exclude_SNs.split(',')
    else:
        exclude_SNs = None

    gps_folder = args.gps_folder
    station_info_file = args.station_info_file
    output_file = args.output_file

    coords = summarize_gps(gps_folder, station_info = station_info_file,
                           output_file = output_file + '.csv', t1 = args.t1, t2 = args.t2,
                           include_SN = args.include_SNs, exclude_SN = args.exclude_SNs)
    response = gemlog.get_gem_response(gain = 'high')
    if 'station' not in coords.keys():
        coords['station'] = coords['SN']
    if 'elevation' not in coords.keys():
        coords['elevation'] = [-9999 for SN in coords.SN]
    if 'network' not in coords.keys():
        coords['network'] = ['' for SN in coords.SN]
    if 'location' not in coords.keys():
        coords['location'] = ['' for SN in coords.SN]
        
    coords['sn'] = coords['SN'] # avoids an error in make_gem_inventory
        
    inv = gemlog.make_gem_inventory(coords, coords, response)
    inv.write(output_file + '.xml', format='STATIONXML') 
    inv.write(output_file + '.kml', format = 'KML') # stations only; locs don't get separate points
    
    return
