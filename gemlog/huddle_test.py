import pandas as pd
import numpy as np
import os, glob, obspy 

def unique(list1):
    unique, index = np.unique(list1, return_index=True)
    return sorted(unique)

def verify_huddle_test(path):

    ##path = '.'
    failures = []
    ## Huddle test performance requirements:
    #### >=3 loggers must have barbs facing each other and all within 15 cm, in a turbulence-suppressed semi-enclosed space, sitting on a shared hard surface on top of padding, with good GPS signal, in a site that is not next to a continuous noise source (duty cycle < 80%). Loggers should all start and stop acquisition within 1 minute of each other, and run for at least one week.

    metadata_path = path + '/metadata'
    gps_path = path + '/gps'
    mseed_path = path + '/mseed'
    for test_path in [metadata_path, gps_path, mseed_path]:
        try:
            fn = os.listdir(test_path)
        except:
            raise(Exception(test_path + ' not found'))
        
    SN_list = unique([filename[-14:-11] for filename in glob.glob(path + '/gps/*')])
    print('Identified serial numbers ' + str(SN_list))

    
    ## Individual Metadata:
    for SN in SN_list:
        print('Checking metadata for ' + SN)
        metadata = pd.read_csv(path +'/metadata/' + SN + 'metadata_000.txt', sep = ',')

        #### battery voltage must be in reasonable range (1.7 to 15 V)
        if any(metadata.batt > 15) or any(metadata.batt < 1.7):
            failure_message = SN + ': Impossible Battery Voltage'
            print(failure_message)
            failures.append(failure_message)
        else:
            print('Sufficient Battery Voltage')

        #### temperature must be in reasonable range (-20 to 60 C)
        if any(metadata.temp > 60) or any(metadata.temp <-20): #celsius 'Impossible Temperature'
            failure_message = SN + ': Impossible Temperature'
            print(failure_message)
            failures.append(failure_message)
        else:
            print('Sufficient Temperature Range')

        ############################################
        #### A2 and A3 must be 0-3.1, and dV/dt = 0 should be true <1% of record

        if not (all(metadata.A2 >= 0) & all(metadata.A2 <= 3.1)):
            pass # the pass keyword just means do nothing. we need something here for indentation purposes
        
            ##

        if (np.sum(np.diff(metadata.A2) == 0) / (len(metadata.A2) -1 )) > 0.01:
            pass
            ## fail the check
        else:
            pass
            ## pass the check
        ############################################

        if not (all(metadata.A2 >=0) & all(metadata.A2 <= 3.1)):
            failure_message = SN + ': Bad A2'
            print(failure_message)
            failures.append(failure_message) 
        else:
            print('A2 okay')
        #### minFifoFree and maxFifoUsed should always add to 75

        if any((metadata.minFifoFree + metadata.maxFifoUsed) != 75):
           failure_message = SN + ': Impossible FIFO sum'
           print(failure_message)
           failures.append(failure_message)
        else:
            print('FIFO sum is correct')
        
        ############################################
        #### maxFifoUsed should be less than 5 99% of the time, and should never exceed 25
        ############################################

        #### maxOverruns should always be zero 
        if any(metadata.maxOverruns) !=0:
            failure_message = SN + ': Too many overruns!'
            print(failure_message)
            failures.append(failure_message)
        else:
            print('Sufficient overruns')    
        
        #### unusedStack1 and unusedStackIdle should always be above some threshold 
        if any(metadata.unusedStack1 <= 50) or any(metadata.unusedStackIdle <= 30):
            failure_message = SN + ': Inadequate unused Stack'
            print(failure_message)
            failures.append(failure_message)
        else:
            print('Sufficient Stack')

        ############################################
        #### gpsOnFlag should never be on for more than 3 minutes at a time
        ############################################

        #### find time differences among samples with gps off that are > 180 sec
        if any(np.diff(metadata.t[metadata.gpsOnFlag == 0]) > 180): 
            failure_message = SN + ': GPS ran for too long'
            print(failure_message)
            failures.append(failure_message)
        else:
            print('GPS runtime ok')

        ## individual GPS:
        gps = pd.read_csv(path +'/gps/' + SN + 'gps_000.txt', sep = ',')
        lat_deg_to_meters = 40e6 / 360 # conversion factor from degrees latitude to meters
        lon_deg_to_meters = 40e6 / 360 * np.cos(np.median(gps.lat) * np.pi/180) # smaller at the poles

        ############################################
        #### lat and lon should vary by less than 100 m for each logger (could change)
        #### the maximum difference between GPS fix times should never exceed 20 minutes
        ############################################

        ## Individual SN waveform data:
        stream = obspy.read(path +'/mseed/*..' + SN + '..HDF.mseed')
        stream.merge()
        #### trim the stream to exclude the first and last 5 minutes
        #### dp/dt = 0 should occur for <1% of record (e.g. clipping, flatlining)
        #### SKIP FOR NOW: noise spectrum must exceed spec/2
        #### SKIP FOR NOW: 20% quantile spectra should be close to self-noise spec
        #### SKIP FOR NOW: noise spectra of sensors must agree within 3 dB everywhere and within 1 dB for 90% of frequencies

    ## Group metadata:
    #### at every given time, temperature must agree within 2C for all loggers
    #### all loggers' first and last times should agree within 20 minutes

    
    ## Group GPS
    #### all loggers' average lat and lon should agree within 1 m
    #### all loggers' first and last GPS times should agree within 20 minutes

    ## group waveform data data:
    #### length of converted data should match among all loggers
    #### a "coherent window" has all cross-correlation coefficients > 0.9, passes consistency criterion, and has amplitude above noise spec. 90% of coherent windows should have only nonzero lags, and none should have persistently nonzero lags (define).


