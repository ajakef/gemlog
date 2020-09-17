import numpy as np
import pandas as pd
import glob, obspy, os

#response = nrl.get_response(sensor_keys = ['Gem', 'Gem Infrasound Sensor v1.0'], datalogger_keys = ['Gem', 'Gem Infrasound Logger v1.0', '0 - 128000 counts/V'])

def _get_station_info(station_info):
    required_keys = ['SN', 'network', 'station', 'location']
    if type(station_info) == str:
        station_info = pd.read_csv(station_info, names = required_keys, dtype = {key:'str' for key in required_keys})
    elif (type(station_info) is not pd.DataFrame) or any([key not in station_info.keys() for key in required_keys]):
        raise Exception('invalid station_info')
    return station_info

def make_gem_inventory(station_info, coords, response = 'default'):
    ## ensure that the provided station_info has all the necessary info, and fail if it doesn't
    station_info = _get_station_info(station_info)
    ## ensure that the coords info has all the necessary info, and fail if it doesn't
    if (type(coords) is not pd.DataFrame):
        raise Exception('invalid coords; must be pandas.DataFrame')
    elif ('SN' not in coords.keys()) and any([key not in coords.keys() for key in ['network', 'station', 'location']]):
        raise Exception("invalid coords; must contain key 'SN', or keys 'network', 'station', 'location'")
    elif any([key not in coords.keys() for key in ['lat', 'lon']]):
        raise Exception("invalid coords; must contain keys 'lat', 'lon'")
    
    ## create the inventory and loop through all the networks, stations, and locations in it
    inventory = obspy.Inventory()
    all_networks = _unique(station_info['network'])
    for network_name in all_networks:
        network = obspy.core.inventory.Network(network_name)
        ## loop through all the stations in the current network
        all_stations = _unique(station_info[station_info['network'] == network_name]['station'])
        for station_name in all_stations:
            station = obspy.core.inventory.Station(station_name, latitude = 0, longitude = 0, elevation = 0)
            ## loop through all the locations/channels in the current station
            all_locations = _unique(station_info[(station_info['network'] == network_name) & (station_info['station'] == station_name)]['location'])
            for location_name in all_locations:
                ## find the serial number for this location
                SN = station_info[(station_info['network'] == network_name) & (station_info['station'] == station_name) & (station_info['location'] == location_name)]['SN'].iloc[0]
                ## find the coords line corresponding to this location
                if 'SN' in coords.keys():
                    line = coords[coords['SN'] == SN]
                else:
                    line = coords[(coords['network'] == network_name) & (coords['station'] == station_name) & (coords['location'] == location_name)]
                ## We need to extract the coordinate and times for this location.
                lat = line['lat'].iloc[0]
                lon = line['lon'].iloc[0]
                ## Check to see if the coords include start/end times. If so, pad
                ## them 1 hour to be safe. If not, assume the station is eternal.
                if all([key in coords.keys() for key in ['starttime', 'endtime']]):
                    t1 = obspy.UTCDateTime(line['starttime'].iloc[0]) - 3600 
                    t2 = obspy.UTCDateTime(line['endtime'].iloc[0]) + 3600
                else:
                    t1 = obspy.UTCDateTime('1970-01-01T00:00:00')
                    t2 = obspy.UTCDateTime('9999-12-31T23:59:59')
                ## Assume it's a gem v1.0--safe for now
                equipment = obspy.core.inventory.util.Equipment(serial_number = SN, model = 'Gem Infrasound Logger v1.0')
                channel = obspy.core.inventory.Channel('HDF', location_code = location_name, latitude = lat, longitude = lon, elevation = 0, depth = 0, response = response, equipments = equipment, start_date = t1, end_date = t2)
                station.channels.append(channel)
            ## calculate the current station's coordinate as the mean of its locations,
            ## then append it to the network
            station.latitude = np.mean([channel.latitude for channel in station.channels])
            station.longitude = np.mean([channel.longitude for channel in station.channels])
            network.stations.append(station)
        ## append the current network to the inventory
        inventory.networks.append(network)
    return inventory
    
def rename_files(infile_pattern, station_info, output_dir, output_format = 'mseed'):
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
        w = np.where(station_info.SN == SN)[0][0]
        network = station_info.network[w]
        station = station_info.station[w]
        location = station_info.location[w]
        for tr in st:
            tr.stats.network = network
            tr.stats.station = station
            tr.stats.location = location
        outputFile = output_dir + '/' + '%s.%s.%s.%s.HDF.%s' % (fileParts[0], network, station, location, output_format)
        st.write(outputFile, format = output_format)
        print(str(i) + ' of ' + str(len(infiles)) + ': ' + infile + ', ' + outputFile)
    return station_info



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
    w = np.where((np.abs(x.lat - np.median(x.lat)) < (N * np.std(x.lat))) & (np.abs(x.lon - np.median(x.lon)) < (N * np.std(x.lon))))[0]
    if len(w) < x.shape[0]:
        return _remove_outliers(x.iloc[w,:], N=N)
    else:
        return(x)

def read_gps(gps_dir_pattern, SN):
    gpsDirList = sorted(glob.glob(gps_dir_pattern))
    gpsTable = pd.DataFrame(columns=['year', 'date', 'lat', 'lon', 't'])
    for gpsDir in gpsDirList:
        fnList = glob.glob(gpsDir + '/' + SN + '*')
        if len(fnList) > 0: # if any gps files matching SN are found, read and append the last
            gpsTable = gpsTable.append(pd.read_csv(sorted(fnList)[-1]), ignore_index=True)
    return gpsTable
ReadLoggerGPS = read_gps # alias; v1.0.0

def summarize_gps(gps_dir_pattern, output_file = '', station_info = None):
    gpsDirList = sorted(glob.glob(gps_dir_pattern))
    gpsFileList = []
    for gpsDir in gpsDirList:
        gpsFileList += glob.glob(gpsDir + '/' + '*gps*txt')
    snList = []
    for gpsFile in gpsFileList:
        snList.append(gpsFile.split('/')[-1][:3])
    snList = sorted(_unique(snList))
    coords = pd.DataFrame(columns = ['SN', 'lat', 'lon', 'lat_SE', 'lon_SE', 'starttime', 'endtime', 'num_samples'])
    avgFun = lambda x: np.mean(x)
    seFun = lambda x: np.std(x)/np.sqrt(len(x))
    for i, SN in enumerate(snList):
        print(str(i) + ' of ' + str(len(snList)) + ': ' + SN)
        gpsTable = _remove_outliers(read_gps(gps_dir_pattern, SN))
        if gpsTable.shape[0] > 0:
            coords = coords.append(pd.DataFrame(
                [[SN, avgFun(gpsTable.lat), avgFun(gpsTable.lon), seFun(gpsTable.lat), seFun(gpsTable.lon), gpsTable.t.min(), gpsTable.t.max(), gpsTable.shape[0]]],
                columns = ['SN', 'lat', 'lon', 'lat_SE', 'lon_SE', 'starttime', 'endtime', 'num_samples']), ignore_index = True)
    if station_info is not None:
        station_info = _get_station_info(station_info)
        network = []
        station = []
        location = []
        for SN in coords.SN:
            try:
                w = np.where(station_info.SN == SN)[0][0]
                network.append(station_info.network[w])
                station.append(station_info.station[w])
                location.append(station_info.location[w])
            except:
                network.append('')
                station.append('')
                location.append('')
        coords['network'] = network
        coords['station'] = station
        coords['location'] = location
    if output_file is not None:
        coords.to_csv(output_file)
    return coords
    
SummarizeAllGPS = summarize_gps # alias, v1.0.0
