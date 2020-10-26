import pandas as pd
import numpy as np

SN = '077'
metadata = pd.read_csv('metadata/' + SN + 'metadata_000.txt', sep = ',')

## Huddle test performance requirements:
#### >=3 loggers must have barbs facing each other and all within 15 cm, in a turbulence-suppressed semi-enclosed space, sitting on a shared hard surface on top of padding, with good GPS signal, in a site that is not next to a continuous noise source (duty cycle < 80%). Loggers should all start and stop acquisition within 1 minute of each other, and run for at least one week.

## Metadata:
#### battery voltage must be in reasonable range (1.7 to 15 V)
if any(metadata.batt > 15) or any(metadata.batt < 1.7):
    raise(Exception('Impossible Battery Voltage'))
else:
    print('Sufficient Battery Voltage')
#%%
#### temperature must be in reasonable range (-20 to 60 C)
if any(metadata.temp > 60) or any(metadata.temp <-20): #celcius
    raise(Exception('Impossible Temperature'))
else:
    print('Sufficient Temperature Range')
#%%
#### at every given time, temperature must agree within 2C for all loggers
#### A2 and A3 must be 0-3.1, and dV/dt = 0 should be true <1% of record
#%%
#### minFifoFree and maxFifoUsed should always add to 75
if any((metadata.minFifoFree + metadata.maxFifoUsed) != 75):
    raise(Exception('Impossible FIFO sum'))
else:
    print('FIFO sum is correct')
#### maxFifoUsed should be less than 5 99% of the time, and should never exceed 25
#%%
#### maxOverruns should always be zero
if any(metadata.maxOverruns) !=0:
    raise(Exception('Too many overruns!'))
else:
    print('Sufficient overruns')    
#%%
#### unusedStack1 and unusedStackIdle should always be >50 (this could change)
if any(metadata.unusedStack1) or any(metadata.unusedStackIdle) >50:
    raise(Exception('Inadequate unused Stack'))
else:
    print('Sufficient Stack')
#%%
#### gpsOnFlag should never be on for more than 3 minutes at a time
#### find time differences among samples with gps off that are > 180 sec
if any(np.diff(metadata.t[metadata.gpsOnFlag == 0]) > 180):
    raise(Exception('GPS ran for too long'))
else:
    print('GPS runtime ok')
#### all loggers' first and last times should agree within 20 minutes

## GPS:
#### lat and lon should vary by less than 100 m for each logger (could change)
#### all loggers' average lat and lon should agree within 1 m
#### the maximum difference between GPS fix times should never exceed 20 minutes
#### all loggers' first and last GPS times should agree within 20 minutes

## Waveform data (excluding first and last hour)--all of these could change:
#### length of converted data should match among all loggers
#### a "coherent window" has all cross-correlation coefficients > 0.9, passes consistency criterion, and has amplitude above noise spec. 90% of coherent windows should have only nonzero lags, and none should have persistently nonzero lags (define).
#### dp/dt = 0 should occur for <1% of record (e.g. clipping, flatlining)
#### noise spectrum must exceed spec/2
#### noise spectra of sensors must agree within 3 dB everywhere and within 1 dB for 90% of frequencies
#### 20% quantile spectra should be close to self-noise spec (check to see if this holds up)
