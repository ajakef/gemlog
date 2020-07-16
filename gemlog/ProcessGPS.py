import pandas as pd
import numpy as np
import glob, os

# function to get unique values from list while preserving order (set-based shortcut doesn't do this)
def unique(list1): 
    unique_list = []
    for x in list1:
        # check if exists in unique_list or not 
        if x not in unique_list: 
            unique_list.append(x) 
    return unique_list

## function to exclude outliers by repeatedly calculating standard dev and tossing points outside N standard devs, until none are left
## this matters because occasionally a dataset will start or end with short recordings made elsewhere, which must be excluded from the calculation of station coords
def RemoveOutliers(x, N = 5):
    w = np.where((np.abs(x.lat - np.median(x.lat)) < (N * np.std(x.lat))) & (np.abs(x.lon - np.median(x.lon)) < (N * np.std(x.lon))))[0]
    if len(w) < x.shape[0]:
        return RemoveOutliers(x.iloc[w,:], N=N)
    else:
        return(x)

def ReadLoggerGPS(gpsDirPattern, SN):
    gpsDirList = sorted(glob.glob(gpsDirPattern))
    gpsTable = pd.DataFrame(columns=['year', 'date', 'lat', 'lon', 't'])
    for gpsDir in gpsDirList:
        fnList = glob.glob(gpsDir + '/' + SN + '*')
        if len(fnList) > 0: # if any gps files matching SN are found, read and append the last
            gpsTable = gpsTable.append(pd.read_csv(sorted(fnList)[-1]), ignore_index=True)
    return gpsTable


def SummarizeAllGPS(gpsDirPattern, outputFilename = ''):
    gpsDirList = sorted(glob.glob(gpsDirPattern))
    gpsFileList = []
    for gpsDir in gpsDirList:
        gpsFileList += glob.glob(gpsDir + '/' + '*gps*txt')
    snList = []
    for gpsFile in gpsFileList:
        snList.append(gpsFile.split('/')[-1][:3])
    snList = sorted(unique(snList))
    coords = pd.DataFrame(columns = ['SN', 'lat', 'lon', 'lat_SE', 'lon_SE', 't1', 't2', 'count'])
    avgFun = lambda x: np.mean(x)
    seFun = lambda x: np.std(x)/np.sqrt(len(x))
    for i, SN in enumerate(snList):
        print(str(i) + ' of ' + str(len(snList)) + ': ' + SN)
        #gpsTable = pd.DataFrame(columns=['year', 'date', 'lat', 'lon', 't'])
        #for gpsDir in gpsDirList:
        #    fnList = glob.glob(gpsDir + '/' + SN + '*')
        #    if len(fnList) > 0: # if any gps files matching SN are found, read and append the last
        #        gpsTable = gpsTable.append(pd.read_csv(sorted(fnList)[-1]))
        gpsTable = RemoveOutliers(ReadLoggerGPS(gpsDirPattern, SN))
        if gpsTable.shape[0] > 0:
            coords = coords.append(pd.DataFrame(
                [[SN, avgFun(gpsTable.lat), avgFun(gpsTable.lon), seFun(gpsTable.lat), seFun(gpsTable.lon), gpsTable.t.min(), gpsTable.t.max(), gpsTable.shape[0]]],
                columns = ['SN', 'lat', 'lon', 'lat_SE', 'lon_SE', 't1', 't2', 'count']), ignore_index = True)
    if not stationFile is None:
        stationInfo = pd.read_csv(stationFile, names = ['SN', 'network', 'station', 'location'], dtype = {'network':'str', 'SN':'str', 'station':'str', 'location':'str'})
        network = []
        station = []
        location = []
        for SN in coords.SN:
            try:
                w = np.where(stationInfo.SN == SN)[0]
                network.append(stationInfo.network[w])
                station.append(stationInfo.station[w])
                location.append(stationInfo.location[w])
            except:
                network.append(np.nan)
                station.append(np.nan)
                location.append(np.nan)
        coords['network'] = network
        coords['station'] = station
        coords['location'] = location
    if len(outputFilename) > 0:
        coords.to_csv(outputFilename)
    return coords
    
