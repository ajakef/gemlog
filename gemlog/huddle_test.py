import pandas as pd
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os, glob, obspy, gemlog
import time
import datetime
from gemlog.gemlog_aux import check_lags
from io import StringIO 
import sys
import pdb
import shutil
from fpdf import FPDF
from matplotlib.backends.backend_pdf import PdfPages

## TO DO:
# - graphs for GPS runtime
# - output to pdf (https://matplotlib.org/stable/gallery/misc/multipage_pdf.html#sphx-glr-gallery-misc-multipage-pdf-py)
# - one plot for each Sn with three axis (normalize GPS and plot with battery voltage)
# - average voltage decay rate

def unique(list1):
    unique, index = np.unique(list1, return_index=True)
    return sorted(unique)

def verify_huddle_test(path, SN_list = [], SN_to_exclude = [], individual_only = False, run_crosscorrelation_checks = False):
    """Perform a battery of tests on converted data from a huddle test to ensure that no Gems are
    obviously malfunctioning. 

    verify_huddle_test assumes that the Gems start recording data simultaneously, stop recording 
    data simultaneously, and are all located in the same place (barbs within about 30 cm). 
    
    Checks performed on each Gem individually include the following:
    --Reported temperature, battery power, and A2/A3 are reasonable. Unreasonable values could 
    indicate hardware problems.
    --Other metadata show no evidence of memory overruns or other computational problems.
    --GPS information is collected frequently and shows a consistent location.
    --Waveforms do not appear to be clipped or otherwise problematic.

    Checks performed by comparing Gems to each other for consistency include the following:
    --All Gems record approximately the same temperature time series.
    --GPS data show all Gems being in approximately the same place.
    --During high signal/noise ratio periods, time lags among sensors are mostly zero. If this is 
    not true, it means that sample timing is being calculated badly, so array analyses could be
    inaccurate.

    Parameters
    ----------
    path : str
        Path to project folder, which contains subfolders mseed, metadata, and gps
    SN_list : list
        List of serial numbers that should be processed
    SN_to_exclude : list
        List of serial numbers that should NOT be processed

    Returns
    -------
    dict with following elements:
    --errors: list of problems that must be fixed
    --warnings: list of potential problems
    --stats: data frame showing quantitative results for all tests
    --results: data frame showing qualitative results for all tests
    """
#%%  
    if os.getlogin() == 'tamara':
        path = '/home/tamara/gemlog/demo_QC'
    elif os.getlogin() == 'jake':
        path = '/home/jake/Work/gemlog_python/demo_QC'
    else:
        print('unknown user, need to define path')
    SN_list = ['058','061','065','077']
    SN_to_exclude = []
    individual_only = False
    run_crosscorrelation_checks = False
#%%    
    ## Create folder for figures
    figure_path = os.path.join(path, "figures")
    file_exists = os.path.exists(figure_path)
    if file_exists == False:    
        os.mkdir(figure_path)
    else:
        pass
    
    
    ## Create blank error, warning, and notes lists to save console messages
    errors = []
    warnings = []
    notes = []
    
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

    ## identify the list of serial numbers to test, including user input to include/exclude
    if len(SN_list) == 0:
        SN_list = unique([filename[-14:-11] for filename in glob.glob(path + '/gps/*')])
    SN_list = [SN for SN in SN_list if SN not in SN_to_exclude]
    print('Identified serial numbers ' + str(SN_list))
   
    pstats_df = pd.DataFrame(index = SN_list) #create dataframe for parameter statistics
    errors_df = pd.DataFrame(index = SN_list) #create dataframe for errors by category
    group_df = pd.DataFrame(index = SN_list) #create dataframe for group gps test
    
    gps_dict = {} 
    metadata_dict = {}
   
    #### Individual Metadata: 
    ### Initialize plots
    plt.close('all') #close any previous figures  
    
    ## Create a PDF output of plots
    # maybe delete?
    plot_dir = 'plots'
    
            
    ##Create battery and temperature time series graphs for all SN
    batt_temp_fig = plt.figure(0,figsize = (6.5,5))
    batt_temp_ax = batt_temp_fig.subplots(2)
    batt_temp_ax[0].set_title("Battery Voltage")
    batt_temp_ax[0].set_ylabel("voltage (V)")
    batt_temp_ax[1].set_title("Temperature")
    batt_temp_ax[1].set_ylabel("temperature (C)")
    batt_temp_ax[1].set_xlabel("month-date hour")
    batt_temp_fig.tight_layout(h_pad = 3.0)
    
    ##Create A2 and A3 time series graphs for all SN
    A2_A3_fig = plt.figure(1, figsize = (6.5,5))
    A2_A3_ax = A2_A3_fig.subplots(2)
    A2_A3_ax[0].set_title("A2")
    A2_A3_ax[0].set_ylabel("Voltage (V)")
    A2_A3_ax[0].set_xlabel("month-date hour")
    A2_A3_ax[1].set_title("A3")
    A2_A3_ax[1].set_ylabel("Voltage (V)")
    A2_A3_ax[1].set_xlabel("month-date hour")
    A2_A3_fig.tight_layout(h_pad = 3.0)
    
    #Create GPS runtime histogram plots for all SN 
    # will not plot for single SN
    # get list of axis even if there is only one in SN_list
    gps_fig = plt.figure(2, figsize = (6.5,5),)
    gps_ax = gps_fig.subplots(len(SN_list))
    gps_ax[0].set_title("GPS Runtime")
    gps_fig.tight_layout()

    
    
        ## Individual Metadata tests:
    for SN_index, SN in enumerate(SN_list):
        print('\nChecking metadata for ' + SN)
        metadata = pd.read_csv(path +'/metadata/' + SN + 'metadata_000.txt', sep = ',')
        metadata_dict[SN] = metadata
        if len(metadata) < 86400: 
            formatter = mdates.DateFormatter('%H:%M')
            xlabel = np.round(metadata.iloc[1,13])
            xlabel = datetime.datetime.utcfromtimestamp(xlabel).strftime('%Y-%m-%d')
        else:
            formatter = mdates.DateFormatter('%m-%d')
            xlabel = "Month-Date"
        
        interval = np.mean(np.diff(metadata.t)) # calculate interval between metadata sampling
        
        ### Battery voltage must be in reasonable range
        # Define battery voltage range
        batt_min = 1.7
        batt_max = 15
        
        
        # Save a few statistics from battery metadata
        pstats_df.loc[SN, "battery min"] = min(metadata.batt)
        pstats_df.loc[SN,"battery max"] = max(metadata.batt)
        
        # Battery voltage minimum tests
        if pstats_df.loc[SN, "battery min"] < batt_min:
            errors_df.loc[SN, "battery min"] = "ERROR"
            err_message = f"{SN} BATTERY ERROR: battery level {np.round(1.7-min(metadata.batt),decimals=2)} Volts below minimum threshold (1.7V)."
            print(err_message)
            errors.append(err_message)
        elif pstats_df.loc[SN, "battery min"] < batt_min + 1.5:
            errors_df.loc[SN, "battery min"] = "WARNING"
            warn_message = f"{SN} BATTERY WARNING: low battery level within {np.round(min(metadata.batt-1.7),decimals=2)} Volts of minimum threshold (1.7V)."
            print(warn_message)
            warnings.append(warn_message)
        else:
            errors_df.loc[SN, "battery min"] = "OKAY"
            
        # Battery voltage maximum tests    
        if pstats_df.loc[SN, "battery max"] > batt_max:
            errors_df.loc[SN, "battery max"] = "ERROR"
            err_message = f"{SN} BATTERY ERROR: battery level {np.round(max(metadata.batt)-15,decimals=2)} Volts above maximum threshold (15V)."
            print(err_message)
            errors.append(err_message)
        elif pstats_df.loc[SN, "battery max"] > batt_max - 0.05:
            errors_df.loc[SN, "battery max"] = "WARNING"
            warn_message = f"{SN} BATTERY WARNING: battery level within {np.round(15 - max(metadata.batt-1.7), decimals=2)} Volts of maximum threshold (15V)."
            print(warn_message)
            warnings.append(warn_message)
        else:
            errors_df.loc[SN, "battery max"] = "OKAY"  
            
        #slope of voltage decay
        
    ##%%%%%##
        ## Temperature must be within reasonable range
        # Define temperature range
        temp_min = -20
        temp_max = 60
        
        # Save a few statistics from temperature metadata
        pstats_df.loc[SN, "temperature min"] = min(metadata.temp)
        pstats_df.loc[SN, "temperature max"] = max(metadata.temp)
        pstats_df.loc[SN, "temperature average"] = np.mean(metadata.temp)
     
        # Temperature minimum check
        if pstats_df.loc[SN,"temperature min"] < -20: #degrees Celsius
            errors_df.loc[SN, "temperature min"] = "ERROR"
            err_message = f"{SN} TEMPERATURE ERROR: temperature {np.abs(np.round(min(metadata.temp)+20,decimals=2))} degrees below minimum threshold (-20 C)."
            print(err_message)
            errors.append(err_message)
        elif pstats_df.loc[SN,"temperature min"] < -15: #modify as needed, just a backbone structure for now.
            errors_df.loc[SN, "temperature min"] = "WARNING"
            warn_message = f"{SN} TEMPERATURE WARNING: temperature within {np.abs(np.round(20 + min(metadata.temp),decimals=2))} degrees of minimum threshold (-20 C)"
            print(warn_message)
            warnings.append(warn_message)
        else:
            errors_df.loc[SN, "temperature min"] = "OKAY" 
        # Temperature maximum check
        if pstats_df.loc[SN, "temperature max"] > 60: #degrees Celsius
            errors_df.loc[SN, "temperature max"] = "ERROR"
            err_message = f"{SN} TEMPERATURE ERROR: temperature {np.round(max(metadata.temp)-60,decimals=2)} degrees above threshold (60 C)."
            print(err_message)
            errors.append(err_message)
        elif pstats_df.loc[SN, "temperature max"] > 50: #modify as needed, just a backbone structure for now.
            errors_df.loc[SN, "temperature max"] = "WARNING"
            warn_message = f"{SN} TEMPERATURE WARNING: temperature within {np.round(60-max(metadata.temp),decimals=2)} degrees of maximum threshold (60 C)."
            print(warn_message)
            warnings.append(warn_message)
        else:
            errors_df.loc[SN,"temperature max"] = "OKAY" 
            
    ##%%%%%##    
        ###Create plots for battery voltage and temperature
        dec_factor = 10 # decimation factor
        
        ##format data for plotting
        #decimate to speed up plotting time
        batt_ind = np.arange(len(metadata.batt)/dec_factor)*dec_factor
        batt_dec = metadata.batt[batt_ind]
        time_unix = round(metadata.t,0)
        time_unix_ind = np.arange(len(time_unix)/dec_factor)*dec_factor
        time_unix_dec = time_unix[time_unix_ind]
        time_datestamp_dec = [datetime.datetime.utcfromtimestamp(int(t)) for t in time_unix_dec]
        temp_ind = np.arange(len(metadata.temp)/dec_factor)*dec_factor
        temp_dec = metadata.temp[temp_ind]
        
        #Plot battery voltage
        batt_temp_ax[0].plot(time_datestamp_dec, batt_dec, label= int(SN))
        batt_temp_ax[0].legend()
        batt_temp_ax[0].xaxis.set_major_formatter(formatter)
        batt_temp_ax[0].set_xlabel(xlabel)
        batt_temp_ax[0].axhline(batt_min, color="red", linestyle = ":", linewidth = 1)
        batt_temp_ax[0].axhline(batt_max, color="red", linestyle=":", linewidth = 1)

    
        #Plot temperature
        batt_temp_ax[1].plot(time_datestamp_dec, temp_dec)
        batt_temp_ax[1].legend(SN_list)
        batt_temp_ax[1].xaxis.set_major_formatter(formatter)
        batt_temp_ax[1].set_xlabel(xlabel)
        
        batt_temp_fig_path = f"{path}/figures/batt_temp.png"
        batt_temp_fig.savefig(batt_temp_fig_path, dpi=300)
                  
    ##%%%%%##         
        #### A2 and A3 must be within a specified range, and dV/dt = 0 should be true <1% of record
        # Define A2 and A3 range
        A_min = 0
        A_max = 3.1 # [TROUBLESHOOT]: May need to be re-evaluated
         
        ## A2
        # Calculate the proportion where A2 is zero
        A2_zerodiff_proportion = (np.sum(np.diff(metadata.A2) == 0) / (len(metadata.A2) -1 )) 
        pstats_df.loc[SN, "A2 flat"] = A2_zerodiff_proportion 
        # True/False test to see whether all data is within range
        within_A2_range = (all(metadata.A2 >= A_min) & all(metadata.A2 <= A_max))
        pstats_df.loc[SN, "A2 range"] = within_A2_range
        
        #Check that A2 dV/dt == 0 less than 99% of time. >99% indicates a likely short circuit.
        if pstats_df.loc[SN, "A2 flat"] > 0.99: 
            errors_df.loc[SN, "A2 flat"] = "ERROR"
            err_message = f"{SN} A2 ERROR: {np.round(A2_zerodiff_proportion*100,decimals=1)}% of A2 dV/dt is exactly 0. More than 99% indicates a likely short circuit."
            print(err_message)
            errors.append(err_message)
        elif pstats_df.loc[SN,"A2 flat"] > 0.95: #95% placeholder; uncertain interpretation
            errors_df.loc[SN, "A2 flat"] = "WARNING"
            warn_message = f"{SN} A2 WARNING: {np.round(A2_zerodiff_proportion*100,decimals=1)}% of A2 dV/dt is exactly 0. More than 99% indicates a likely short circuit."
            warnings.append(warn_message)
            print(warn_message)
        else:
            errors_df.loc[SN, "A2 flat"] = "OKAY" 
            
       #Check A2 is within range
        if not within_A2_range:
            pstats_df.loc[SN, "A2 range"] = 0 #outside of range (false)
            errors_df.loc[SN, "A2 range"] = "ERROR"
            err_message = f"{SN} A2 ERROR: A2 outside of range"
            errors.append(err_message)
            print(err_message)
        else:
            pstats_df.loc[SN,"A2 range"] = 1 #within range(true)
            errors_df.loc[SN,"A2 range"] = "OKAY" 
    ##%%%%%## 
        ## A3
        # Calculate the proportion where A3 is zero
        A3_nonzero = (np.sum(np.diff(metadata.A3) == 0) / (len(metadata.A3) -1 ))
        pstats_df.loc[SN,"A3 flat"] = A3_nonzero
        #True/False test to see whether all data is within range
        within_A3_range = (all(metadata.A3 >= A_min) & all(metadata.A3 <= A_max))
        pstats_df.loc[SN, "A3 range"] = within_A3_range
        
        #Check that A3 dV/dt == 0 less than 99% of time
        if pstats_df.loc[SN,"A3 flat"] > 0.99: 
            errors_df.loc[SN,"A3 flat"] = "ERROR"
            err_message = f"{SN} A3 ERROR: {np.round(A3_nonzero*100,decimals=1)}% of A3 dV/dt is exactly 0. More than 99% indicates a likely short circuit."
            errors.append(err_message)
            print(err_message)
        elif pstats_df.loc[SN,"A3 flat"] > 0.95: #placeholder of 95%
            errors_df.loc[SN, "A3 flat"] = "WARNING"
            warn_message = f"{SN} A3 WARNING: {np.round(A3_nonzero*100,decimals=1)}% of A3 dV/dt is exactly 0. More than 99% indicates a likely short circuit."
            warnings.append(warn_message)
            print(warn_message)
        else:
            errors_df.loc[SN,"A3 flat"] = "OKAY"
            
       #Check A3 is within range     
        if not within_A3_range:
            pstats_df.loc[SN, "A3 range"] = 0 #outside of range (false)
            errors_df.loc[SN, "A3 range"] = "ERROR"
            err_message = f"{SN} A3 ERROR: A3 outside of range."
            errors.append(err_message)
            print(err_message)
        else:
            pstats_df.loc[SN, "A3 range"] = 1 #within range(true)
            errors_df.loc[SN, "A3 range"] = "OKAY" 
            
    ##%%%%%## 
        ##format data to plot A2 and A3 metadata
        A2_ind = np.arange(len(metadata.A2)/dec_factor)*dec_factor
        A2_dec = metadata.A2[A2_ind]
        A3_ind = np.arange(len(metadata.A3)/dec_factor)*dec_factor
        A3_dec = metadata.A3[A3_ind]
        
        #Plot A2 data
        A2_A3_ax[0].plot(time_datestamp_dec, A2_dec)  
        A2_A3_ax[0].legend(SN_list)
        A2_A3_ax[0].xaxis.set_major_formatter(formatter)
        A2_A3_ax[0].set_xlabel(xlabel)
        #Plot A3 data
        A2_A3_ax[1].plot(time_datestamp_dec, A3_dec)
        A2_A3_ax[1].legend(SN_list)
        A2_A3_ax[1].xaxis.set_major_formatter(formatter)
        A2_A3_ax[1].set_xlabel(xlabel)
        
        A2_A3_fig_path = f"{path}/figures/A2_A3.png"
        A2_A3_fig.savefig(A2_A3_fig_path, dpi = 300)
        
    ##%%%%%##            
        #### minFifoFree and maxFifoUsed should always add to 75
        # Determine threshold
        fifo_rule = 75
        fifo_sum = metadata.minFifoFree + metadata.maxFifoUsed
        pstats_df.loc[SN, "fifo sum"] = max(fifo_sum)
        if any((fifo_sum) != fifo_rule):
           errors_df.loc[SN, "fifo sum"] = "INFO"
           notes_message = f"{SN} fifo notes: fifo sum exceeds range by {np.round(fifo_rule - max(fifo_sum),decimals=2)}."
           notes.append(notes_message)
           #print(notes_message)
        ##elif for warning??
        else:
            errors_df.loc[SN,"fifo within range"] = "OKAY"
            
        #### maxFifoUsed should be less than 5 99% of the time, and should never exceed 25
        #how to display in dataframe if less than 5 99%?
        max_fifo_check = np.sum(metadata.maxFifoUsed > 5)/(len(metadata.maxFifoUsed) -1 )
        pstats_df.loc[SN,"max fifo within range"] = max_fifo_check
        if max_fifo_check > 0.01:
            errors_df.loc[SN,"max fifo within range"] = "INFO"
            notes_message = f"{SN} max fifo info: max fifo (max_fifo_check - 0.01)*100 % outside the standard operating range."
            #print(notes_message)
            notes.append(notes_message)
        else:
            errors_df.loc[SN, "max fifo within range"] = "OKAY"
            
        if any(metadata.maxFifoUsed > 25):
            errors_df.loc[SN, "max fifo"] = "NOTE"
            notes_message = f"{SN} max fifo note: max fifo exceeds maximum by {np.round(max(metadata.maxFifoUsed)-25,decimals=2)}."
            notes.append(notes_message)
        else:
            errors_df.loc[SN, "max fifo within range"] = "OKAY"
            
    ##%%%%%##             
        #### maxOverruns should always be zero 
        pstats_df.loc[SN, "max overruns"] = max(metadata.maxOverruns)
        if any(metadata.maxOverruns) !=0:
            errors_df.loc[SN, "max overruns"] = "NOTE"
            notes_message = f"{SN} overruns note: maximum overruns does equal 0. ({max(metadata.maxOverruns)})"
            print(notes_message)
            notes.append(notes_message)
        else:
            errors_df.loc[SN, "max overruns"] = "OKAY"
            
    ##%%%%%##             
        #### unusedStack1 and unusedStackIdle should always be above some threshold 
        #pstats_df.loc[SN,"unused stack1 max"] = max(metadata.unusedStack1)
        #pstats_df.loc[SN,"unused stack idle max"] = max(metadata.unusedStackIdle)
        if any(metadata.unusedStack1 <= 30) or any(metadata.unusedStackIdle <= 30):
            errors_df.loc[SN, "unused stack"] = "WARNING"
            warn_message = f"UNUSED STACK WARNING: One value of unused stack exceeds maximum by {np.round(max(metadata.unusedStack1)-30,decimals=2)}."
            print(warn_message)
            warnings.append(warn_message)
        else:
            errors_df.loc[SN, "unused stack"] = "OKAY"

    ##%%%%%##             
        #### find time differences among samples with gps off that are > 180 sec
        # subtract interval from gps_time_check
        gps_time_check = np.diff(metadata.t[metadata.gpsOnFlag == 0])[10:] # skip the first ten seconds and first GPS cycle, which is very long by design
        gps_mean = gps_time_check[gps_time_check > 11] 
        pstats_df.loc[SN, "mean gps run time"] = np.mean(gps_mean)
        if any(gps_time_check > 180): 
            errors_df.loc[SN, "gps run time"] = "WARNING"
            warn_message = f"{SN} GPS WARNING: GPS runtime is {max(gps_time_check)}."
            print(warn_message)
            warnings.append(warn_message)
            
        ## individual GPS:
            #plot GPS histogram for runtime
            #change plotting options 
        gps_on = metadata.gpsOnFlag
        time_cycle = 900 # seconds in cycle (15 minutes)
        gps_proportion = np.sum(gps_on)/len(gps_on) * time_cycle # time proportion that GPS is on
        if gps_proportion > 180:
            gps_proportion = 180
        pstats_df.loc[SN, "on time proportion"] = gps_proportion
        #IDEAS - change to stacked bar chart with different colors for serial numbers
        # how to select colors automatically?
        time_filt = gps_time_check[gps_time_check > 11] - interval
        time_filt[time_filt > 180] = 180
        binsize = np.arange(10,180,10)
        bins = gps_ax[SN_index].hist(time_filt, bins=np.arange(10,180,5))
        y_scale = np.round((max(bins[0])/2) + 1)
        gps_ax[SN_index].errorbar(gps_proportion, y_scale, yerr= y_scale, ecolor = 'r')
        gps_ax[SN_index].set_ylabel('#' + SN_list[SN_index])
        gps_ax[SN_index].set_xticks([]) # not working
        gps_ax[SN_index].axes.xaxis.set_ticklabels([])
        gps_ax[SN_index].xaxis.set_major_locator(plt.MultipleLocator(20))
        gps_ax[SN_index].xaxis.set_minor_locator(plt.MultipleLocator(10))
        if SN == SN_list[-1]:    
            gps_ax[SN_index].set_xlabel('seconds')
            gps_ax[SN_index].axes.xaxis.set_ticklabels([0,20,40,60,80,100,120,140,160,'>180'])
            gps_ax[SN_index].annotate('on time proportion', (gps_proportion, 10), xytext = (gps_proportion + 10 , 15), color = 'r',
                                  arrowprops = dict(arrowstyle = '->', connectionstyle = "angle, angleA = 90, angleB = 0, rad = 10", color = 'r'))
        gps_fig_path = f"{path}/figures/gps_runtime.png"
        gps_fig.savefig(gps_fig_path, dpi=300)
        
           
        gps = pd.read_csv(path +'/gps/' + SN + 'gps_000.txt', sep = ',')
        gps.t = gps.t.apply(obspy.UTCDateTime)
        gps_dict[SN] = gps
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
    
    #Do not omit rows and columns when displaying in console
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', None)
    print('Results')
    print(errors_df)
    print('Stats:')
    print(pstats_df)
    
    
    print("\nSerial number tests complete.") 
    
# ============================================================================= 
    ## Group metadata tests:    
    ## If we're at this point, we have data from multiple loggers and are clear 
    # to run group tests.
# =============================================================================
    ## Before running the group tests, ensure that we actually have data more than one Gem!
    ## If not, add a warning, and return without conducting group tests.
    group_err = []
    group_warnings = []
    group_notes = []
    
    print("\n\nRunning group GPS tests")
    if (len(SN_list) == 1) or individual_only:
        warn_message = 'Test only includes one logger; cannot run comparison tests'
        group_err.append(warn_message)
        #return {'errors':errors, 'warnings':warnings, 'notes':notes}
        
#### all loggers' first and last times should agree within 20 minutes
    max_mins = 20
    max_secs = max_mins * 60
    # Record all start and stop times in group test dataframe
    for SN_index, SN in enumerate(SN_list):
        start_time = min(metadata.t)
        stop_time = max(metadata.t)
        group_df.loc[SN, "start time"] = start_time
        group_df.loc[SN, "end time"] = stop_time
    
    # Find the median and create a range within the max time distance for start and stop times
    gps_start_med = np.median(group_df.iloc[:,0]) # find the median of start times
    gps_stop_med = np.median(group_df.iloc[:,1]) # find the median of stop times
    upper_start = gps_start_med + max_secs/2 # create an upper bound for start times
    lower_start = gps_start_med - max_secs/2 # create a lower bound for start times
    upper_stop = gps_stop_med + max_secs/2 # create an upper bound for stop times
    lower_stop = gps_stop_med - max_secs/2 # create a lower bound for stop times
    
    # Check all SN start and stop times to ensure they are within range
    for SN_index, SN in enumerate(SN_list):
        if not lower_start <= group_df.iloc[SN_index,0] <= upper_start:
            err_message = (f"{SN} GROUP GPS ERROR: The start times are not within {max_mins} minutes of the median.")
            group_err.append(err_message)
            print(err_message)
        if not lower_stop <= group_df.iloc[SN_index,1] < upper_stop:
            err_message = (f"{SN} GROUP GPS ERROR: The stop times are not within {max_mins} minutes of the median.")
            group_err.append(err_message)
            print(err_message)
    if len(group_err) == 0:
        note = "The start and stop times for all loggers agree."
        group_notes.append(note)
        print(note)

#### at every given time, temperature must agree within 2C for all loggers
    interval = 60 #seconds (time between checks)
    
    # Determine start and stop times (Unix time)
    mod = upper_start % interval # modulus remainder to round start time to even minute
    temp_start = upper_start - mod # start time on an even minute
    mod = upper_stop % interval
    temp_end = upper_stop - mod # stop time on an even minute
    
    # Create dataframe to house temperatures to check
    column_index = np.arange(0,int((temp_end-temp_start)/interval),1) # theoretically, the number of values between start and end
    
    group_temp_df = pd.DataFrame(index = SN_list, columns = column_index ) # create dataframe to house temperatures at each minute for each SN
    diff = [] # contain difference calculations
    times_checked = np.zeros((1,max(column_index)))
    
    for SN_index, SN in enumerate(SN_list):
        metadata = metadata_dict[SN]
        argstart_list = [] # reset closest start list
        argend_list = [] # reset closest end list
        times = metadata.t # call all the time metadata
        for time in times: # for each time value in time metadata
            argstart_list.append(np.abs(temp_start - time)) # create a list of time values minus start time (find closest to start)
            argend_list.append(np.abs(temp_end - time)) # create list of time values minus end time (find closest index to end)
        
        start_index = np.argmin(argstart_list) # find index for closest minute time for start
        #format to not include nans from metadata
        stop_index = np.argmin(argend_list) # find index for last time value
        

        temp_times = np.arange(temp_start, temp_end, 60) # for label
        
        times_to_check_index = np.arange(start_index, stop_index, 60) # create evenly spaced array of even minutes
        # must create index based on mutally agreed start time
        for df_index, index in enumerate(times_to_check_index):
            #TROUBLESHOOT: 061 and 065 saving into dataframe as nan after index 29
            group_temp_df.iloc[SN_index,df_index] = metadata.temp[index] # save minute temperature reading into dataframe
            times_checked[0, df_index] = (metadata.t[index]) # ***not efficient***
    error = False
    for column in column_index: # column represents temperature data for each time that will be checked 
        # might be operating dataframe functions on entire set, not by columns...
        temp_median = np.round(group_temp_df[column].median(),2)
        temp_range = np.round(group_temp_df[column].max() - group_temp_df[column].min(),2)
        outliers = (np.where(any(group_temp_df[column]) > temp_median + 1 or any(group_temp_df[column] < temp_median - 1))[0])
        # Find offending serial numbers outside of temperature range
        x = (group_temp_df.index[group_temp_df[column] > temp_median + 1].tolist())       
        if len(x) > 0:
            error = True
            ts = int(times_checked[0,column])
            time_lookup = datetime.datetime.utcfromtimestamp(ts)
            err_message = (f"SN {x} recorded temperatures greater than 1 on either side of the temperature median {temp_median} on {time_lookup}. Total temperature range = {temp_range}")
            group_err.append(err_message)
            print(err_message)
    if error == False:
        print("The recorded temperatures are within two degrees Celcius")  
     
            
    #  Relict temperature check (KeyError: 'SN')      
    # all_temperatures = np.zeros((len(SN_list), len(times_to_check_index)))
    # temperatures = {}
    # for i, SN in enumerate(SN_list):
    #     temperatures[SN] = scipy.signal.lfilter(np.ones(100)/100, [1],
    #         scipy.interpolate.interp1d(metadata_dict[SN].t, metadata_dict[SN].temp)(times_to_check_index))
    #     # to do: use the median instead
    #     # average_temperatures += temperatures[SN]/ len(SN_list)
    #     all_temperatures[i,:] =  temperatures[SN]
    # med_temperatures = np.median(all_temperatures, 0)
    # for SN in SN_list:
    #     if np.sum(np.abs(temperatures[SN] - med_temperatures) > 2)/len(times_to_check_index) > 0.1:
    #         failure_message = SN + ': disagrees excessively with average temperature'
    #         print(failure_message)
    #         errors.append(failure_message)
    #     else:
    #         print(SN + ': Temperatures agree')
   
    #### all loggers' average lat and lon should agree within 1 m

    ## group waveform data data:
    #### length of converted data should match among all loggers
    #### a "coherent window" has all cross-correlation coefficients > 0.9, passes consistency criterion, and has amplitude above noise spec. 90% of coherent windows should have only nonzero lags, and none should have persistently nonzero lags (define).
    if run_crosscorrelation_checks:
        DB = gemlog.make_db(path + '/mseed', '*', 'tmp_db.csv')
        DB = DB.loc[DB.station.isin(SN_list),:]
        [t, lag, xc_coef, consistency] = check_lags(DB, winlength=100)
        coherent_windows = (consistency == 0) & (np.median(xc_coef, 0) > 0.8)
        zero_lags = lag[:,coherent_windows]==0
        num_coherent_windows = np.sum(coherent_windows)
        if (np.sum(np.all(zero_lags, 0)) / num_coherent_windows) < 0.8:
            failure_message = 'Time lags are excessively nonzero for coherent time windows'
            print(failure_message)
            group_err.append(failure_message)
        else:
            note = 'Time lags for coherent time windows are mostly/all zero'
            group_notes.append(note)
            print(note)
        
        ## add a figure showing time lags plotted over time
        time_lags_fig = plt.figure(3)
        for i in range(len(lag)):
            plt.plot(np.array(t)-t[0], lag[i])
        plt.legend() ## not sure how to do this
        plt.title('Time lags from cross-correlation')
        plt.ylabel('Time lag (samples)')
        time_lags_fig_path = f"{path}/figures/time_lags.png"
        time_lags_fig.savefig(time_lags_fig_path, dpi = 300)

    #return {'errors':errors, 'warnings':warnings, 'notes':notes}
    
#%%

# ============================================================================
    # Information on tests including thresholds to be printed in the report
# ============================================================================
    info = []
    trouble = []

    ## Battery test information and troubleshooting ##
    info.append(f"""BATTERY TEST INFO: This test is designed to ensure the voltage of each gem is within [{batt_min} to {batt_max} 
                Volts at each recorded instance. Gems within this specified range, but within a range 0.5 Volts from the threshold 
                will result in a warning message. This message is printed into the console and into the automatically generated pdf
                within the metadata folder. A plot is generated to visualize which serial numbers are malfunctioning
                and where they are operating. This plot is also saved in the metadata folder under an automatically generated folder
                called figures.""")
    trouble.append("""BATTERY TROUBLESHOOTING: If you received a battery warning, it is likely the gem was not able to record any 
                   waveform data. This can usually be fixed by changing the batteries. If you received a battery error [INSERT PROBLEM]
        [INSERT TROUBLESHOOT]
        """)
    ## Temperature test information and troubleshooting ##
    info.append(f"""TEMPERATURE TEST INFO: This test ensures the gemlogger is recording ambient air temperatures within an reasonable 
                range. This range is set at {temp_min} to {temp_max} Celsius or {(temp_min * 9/5) + 32} to {(temp_max * 9/5) + 32} 
                Fahrenheit. At temperatures outside this range, the electrical components of the gem could malfunction.
        """)
    trouble.append("""TEMPERATURE TROUBLESHOOTING: If you received a temperature warning, you are approaching the limit of the 
                   temperature range operation for the gemlogger ({temp_min} to {temp_max} C). If this value does not reflect an 
                   accurate ambient air temperature, [INSERT TROUBLESHOOTING]
        """)
#%%
   
# ============================================================================= 
    # Create a PDF output of plots with the date of report, errors warning and
    # notes list, and metadata summary dataframes
    # Create seperate package to reduce gemlog dependencies for detailed report
# =============================================================================
    report_path = os.path.join(path, "reports")
    file_exists = os.path.exists(report_path)
    if file_exists == False:    
        os.mkdir(report_path)
    else:
        pass
## Set up report pages and headings    
    report_date = datetime.datetime.today()
    report_date = report_date.strftime("%Y-%m-%d")
    filename = str("Huddle_test_output_" + report_date)
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font('helvetica', 'B', size=12)
    pdf.cell(0,10, f"Huddle Test Results", border=0, ln=0, align= 'C')
    pdf.ln() #new line
    pdf.cell(0,10, f"Date: {report_date}",border=0,align= 'C', ln=1)
    pdf.ln()

## Insert error and warning list
    pdf.cell(0,5,"Errors and Warnings", align = 'L', ln=1)
    pdf.set_font('helvetica', size=8)
    for i, error in enumerate(errors):
        pdf.cell(12,4, '%s' % errors[i], ln=1)
    for j, warning in enumerate(warnings):
        pdf.cell(12,4, '%s' % warnings[j], ln=1)
    pdf.ln()    
## Insert notes list       
    pdf.cell(0,5,"Notes", align = 'L', ln=1)
    pdf.set_font('helvetica', size=8)
    for k, note in enumerate(notes):
        pdf.cell(12,3, '%s' % notes[i], ln=1)
    pdf.ln()
    
## Insert errors dataframe as status
    #Format for multiline headers
    status_header = []
    for index, col in enumerate(errors_df.columns):
       col_header = [[],[]]
       col_header = col.split() 
       status_header.append(col_header)
    #reduce excessive word sizes
    for word in status_header:
        if word[0] == "temperature":
            word[0] = "temp"
        if word[1] == "overruns":
            word[1] = "overrun"
            
    space = ''
    for i, header in enumerate(status_header):
        if len(status_header[i]) == 3:
            status_header[i].insert(0,space)
        elif len(status_header[i]) == 2:
            status_header[i].insert(0,space)
            status_header[i].insert(1,space)
        elif len(status_header[i]) == 1:
            status_header[i].insert(0,space)
            status_header[i].insert(1,space)
            status_header[i].insert(2,space)
        else:
            pass
    
    
    # Insert title for summary table
    pdf.set_font('helvetica', 'BU', size=10)
    pdf.cell(0,10, 'Gem Sensor Status Summary', align='C')
    pdf.ln()
    
    # Insert headers for summary table
    pdf.set_font('helvetica', 'B', size=8)
    for j in range(0,4):
        if j == 2:
            pdf.cell(8,10, 'SN') # Create SN column heading
        else:
            pdf.cell(8,10, ' ') #Other lines are blank in SN column heading
        for i, header in enumerate(status_header):
            pdf.cell(12,5, '%s' % status_header[i][j])
        pdf.ln() 
        
    # print errors data frame into PDF as "status"
    pdf.set_font('helvetica', 'B', size=6)
    for i in range(0,len(errors_df)):
       pdf.cell(8,10, '%s' % SN_list[i])
       for j in range(0,len(errors_df.columns)): 
           pdf.cell(12,10, '%s' % errors_df.iloc[i,j], 1, 0, 'C')
       pdf.ln()
       
    
## Insert pstats_df as details        
    #Format for multiline headers
    pdf_header = []
    for index, col in enumerate(pstats_df.columns):
       col_header = [[],[]]
       col_header = col.split() 
       pdf_header.append(col_header)
    #reduce excessive word sizes
    for word in pdf_header:
        if word[0] == "temperature":
            word[0] = "temp"
        if word[1] == "overruns":
            word[1] = "overrun"
            
    space = ''
    for i, header in enumerate(pdf_header):
        if len(pdf_header[i]) == 3:
            pdf_header[i].insert(0,space)
        elif len(pdf_header[i]) == 2:
            pdf_header[i].insert(0,space)
            pdf_header[i].insert(1,space)
        elif len(pdf_header[i]) == 1:
            pdf_header[i].insert(0,space)
            pdf_header[i].insert(1,space)
            pdf_header[i].insert(2,space)
        else:
            pass
    # Insert title for details table
    pdf.ln()
    pdf.set_font('helvetica', 'BU', size=10)
    pdf.cell(0,10, 'Gem Sensor Metadata Details', align='C')
    pdf.ln()
    
    pdf.set_font('helvetica', 'B', size=8)
    
    
    # Create headers for details table
    for j in range(0,4):
        if j == 2:
            pdf.cell(8,10, 'SN') # Create SN column heading
        else:
            pdf.cell(8,10, ' ') #Other lines are blank in SN column heading
        for i, header in enumerate(pdf_header):
            pdf.cell(12,5, '%s' % pdf_header[i][j])
        pdf.ln() 
    
    #print out pstats_df into pdf
    pdf.set_font('helvetica', size = 8)
    for i in range(0,len(pstats_df)):
       pdf.cell(8,10, '%s' % SN_list[i])
       for j in range(0,len(pstats_df.columns)): 
           pdf.cell(12,10, '%s' % np.round((pstats_df.iloc[i,j]),3), 1, 0, 'C')
       pdf.ln()
       
    
    
    ## Add figures into report
    img_height = 120
    img_width = 176
    
    pdf.image(batt_temp_fig_path, w = img_width , h = img_height) 
    pdf.ln()
    pdf.image(A2_A3_fig_path, w = img_width, h = img_height)
    pdf.ln()
    pdf.image(gps_fig_path, w = img_width, h = 135)
    pdf.ln()
    pdf.ln()
    
## Group test results
    pdf.set_font('helvetica', 'B', size = 10)
    pdf.cell(0,10, f"Group Test Results", border=0, ln=0, align= 'C')
    pdf.ln()
    pdf.set_font('helvetica', size=8)
    for i, error in enumerate(group_err):
        pdf.cell(12,4, '%s' % group_err[i], ln=1)
    for i, note in enumerate(group_notes):
         pdf.cell(12,4, '%s' % group_notes[i], ln=1)
    #return {'errors':errors, 'warnings':warnings, 'stats':pstats_df, 'results':errors_df}
    pdf.ln()
    for i, note in enumerate(info):
        pdf.multi_cell(200,5, '%s' %info[i])
        pdf.ln()
    
    ## add the time lags if they were actually calculated
    if run_crosscorrelation_checks:
        pdf.ln()
        pdf.image(time_lags_fig_path, w = img_width, h = img_height)
        
## Close and name file    
    pdf.output(f"{report_path}/{filename}.pdf")
    

