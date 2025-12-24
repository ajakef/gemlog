#from gemlog.core import _make_empty_gps, _make_empty_metadata, _read_single, _calculate_drift, _make_empty_header, _robust_regress
import numpy as np
import matplotlib.pyplot as plt
import obspy, gemlog
from gemlog.exceptions import CorruptRawFile

def find_gps_discontinuity(fn, min_error = 1/900):
    ## code from _read_several
    #D = np.ndarray([0,2]) # expected number 7.2e5
    #fnList = [fn]
    startMillis = 0
    #header = gemlog.core._make_empty_header(fnList)
    #G = gemlog.core._make_empty_gps()
    #M = gemlog.core._make_empty_metadata()
    L = gemlog.core._read_single(fn, startMillis, require_gps = True)
    #if(L['data'][0,0] < startMillis):
    #    L['metadata'].millis += 2**13
    #    L['gps'].msPPS += 2**13
    #    L['data'][:,0] += 2**13
    
    w = np.where(np.abs(np.diff(L['gps'].t) / np.diff(L['gps'].msPPS)/0.001024 -1) > min_error)[0] # unitless time deviation from 1 sec / 1024 gem_millis
    return(L['gps'].iloc[w,:])

def gps_debug_plot(fn, degree = 3, MAD = 0.01, z = 4, version = '1.10'):
        ## plot a file's deviation from linear GPS drift
    ## code from _read_several
    
    D = np.ndarray([0,2]) # expected number 7.2e5
    fnList = [fn]
    startMillis = 0
    header = gemlog.core._make_empty_header(fnList)
    G = gemlog.core._make_empty_gps()
    M = gemlog.core._make_empty_metadata()
    L = gemlog.core._read_single(fn, startMillis, require_gps = True, version = version)
    if(L['data'][0,0] < startMillis):
        L['metadata'].millis += 2**13
        L['gps'].msPPS += 2**13
        L['data'][:,0] += 2**13
    ## if the millis series is discontinuous, reject the file
    dMillis = np.diff(L['data'][:,0])
    if any(dMillis < 0) or any(dMillis > 1000):
        raise CorruptRawFile(f'{fn} sample times are discontinuous, skipping this file')

    ## code from _calculate_drift
    x = L['gps'].msPPS
    y = L['gps'].t

    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(x, y - 0.001024 * x, 'ro')
    try:
        reg, num_gps_nonoutliers, MAD_nonoutliers, resid, xx, yy = gemlog.core._robust_regress(x, y, degree = degree, MAD = MAD, z = z)
        plt.plot(xx, yy - 0.001024 * xx, 'b.')
        plt.subplot(2,1,2)
        plt.plot(xx, resid, 'b.')
    except:
        xx = None; yy = None; resid = None
    return x, y, xx, yy, resid

def time_adjust(fn_old, fn_new, t1, t2, shift, verbose = False):
    t1 = obspy.UTCDateTime(t1)
    t2 = obspy.UTCDateTime(t2)

    # find newline mode
    with open(fn_old, 'rb') as old:
        for i, line in enumerate(old):
            if line[-2] == 13: # \r character
                newline = '\r\n'
            else:
                newline = '\n'
            break
                
    # open the old and new files
    with open(fn_old, 'r', errors = 'ignore') as old:
        with open(fn_new, 'w', newline = newline) as new:

            # loop through lines in old file
            for i, line in enumerate(old):
                if verbose:
                    #breakpoint()
                    print(line)
                # if it is a GPS line, read the time
                if line[0] == 'G':
                    splitline = line.split(',')
                    t = obspy.UTCDateTime(','.join(splitline[3:9]))
                    # if it's between t1 and t2, adjust the time by shift, and change the line accordingly
                    if (t >= t1) and (t <= t2):
                        t += shift
                        fixed_line = ','.join(splitline[:3] + \
                                              ['%d,%d,%d,%d,%d,%.1f' % \
                                              (t.year, t.month, t.day, t.hour, t.minute, t.second)] + \
                                              splitline[9:])
                        line = fixed_line
                # write the GPS line to the new file
                new.write(line)
    
    # done with file loop. run gps_debug on the old and new to confirm that it's good now
    plt.subplot(2,1,1)
    gps_debug(fn_old)
    plt.subplot(2,1,2)
    gps_debug(fn_new)
    plt.tight_layout()

## older code, doesn't show effects of outlier-dropping well

def gps_debug(fn):
    ## plot a file's deviation from linear GPS drift
    ## code from _read_several
    D = np.ndarray([0,2]) # expected number 7.2e5
    fnList = [fn]
    startMillis = 0
    header = _make_empty_header(fnList)
    G = _make_empty_gps()
    M = _make_empty_metadata()
    L = _read_single(fn, startMillis, require_gps = True)
    if(L['data'][0,0] < startMillis):
        L['metadata'].millis += 2**13
        L['gps'].msPPS += 2**13
        L['data'][:,0] += 2**13
    ## if the millis series is discontinuous, reject the file
    dMillis = np.diff(L['data'][:,0])
    if any(dMillis < 0) or any(dMillis > 1000):
        raise CorruptRawFile(f'{fn} sample times are discontinuous, skipping this file')

    ## code from _calculate_drift

    reg, num_gps_nonoutliers, MAD_nonoutliers, resid, xx, yy = _robust_regress(L['gps'].msPPS, L['gps'].t)
    
    plt.plot(xx/1000, yy - reg[3] - xx * 0.001024, 'k.')
    #plt.plot(xx, xx * (reg[2] - 0.001024) + xx**2 * reg[1] + xx**3 * reg[0])
    plt.xlabel('Seconds of file')
    plt.ylabel('GPS Time Deviation (s)')
    plt.title(f'GPS times: file {fn}')
    return reg, num_gps_nonoutliers, MAD_nonoutliers, resid, xx, yy




## 080
#for num in range(420, 424):
#    plt.figure()
#    time_adjust(f'raw/FILE0{num}.080', f'raw_new/FILE0{num}.080', '2023-08-15T19:22:35', '2023-08-16T00:45:35', -2)

## 088
#for num in range(421, 424):
#    plt.figure()
#    time_adjust(f'raw/FILE0{num}.088', f'raw_new/FILE0{num}.088', '2023-08-15T22:58:37', '2023-08-30T00:45:35', 20)



##############################
## code to look for GPS steps

from scipy.ndimage import median_filter

def check_raw_file_for_gps_steps(fn):
    ## returns True if GPS step detected
    G = gemlog.core._read_single(fn, 0, require_gps = True)['gps']
    x = np.array(G['msPPS'])
    y = median_filter(np.array(G['t']) - x * 0.001024, 3) # medfilt *diameter*
    dy = np.diff(y)
    return np.abs(np.sum(dy[(dy > 0.1) | (dy < -0.1)])) > 0.1 # in case a step is reversed quickly, which happened in FILE0004.356 in the roof GPS test

## test:
#print(check_raw_file_for_gps_steps('/home/jake/Work/gemlog_python/data/test_data/gps_time_discontinuity/FILE0061.200')) # has GPS step; should return True
#print(check_raw_file_for_gps_steps('/home/jake/Work/gemlog_python/data/v1.10/FILE0002.232')) # no GPS step; should return False
#print(check_raw_file_for_gps_steps('/home/jake/2024-04-26_RoofTestGPS/raw/FILE0004.356')) # has a reversed GPS step; should return False


## run on GPS roof test data
#import glob
#filenames = sorted(glob.glob('*'))
#step_files = []
#for fn in filenames:
#    has_steps = check_raw_file_for_gps_steps(fn)
#    if has_steps:
#        step_files.append(fn)
#    print((fn, check_raw_file_for_gps_steps(fn)))


###################################



    
