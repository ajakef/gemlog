import pandas as pd
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import os, glob, obspy, gemlog
import time
import datetime
from gemlog.gemlog_aux import check_lags
from io import StringIO 
import sys
import pdb
import shutil
from fpdf import FPDF

## TO DO:
# - graphs for GPS runtime
# - output to pdf (https://matplotlib.org/stable/gallery/misc/multipage_pdf.html#sphx-glr-gallery-misc-multipage-pdf-py)
# - one plot for each Sn with three axis (normalize GPS and plot with battery voltage)

def unique(list1):
    unique, index = np.unique(list1, return_index=True)
    return sorted(unique)

def verify_huddle_test(path, SN_list = [], SN_to_exclude = [], individual_only = False):
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
    
    # SN_list = ['058','061','065','077']
    # SN_to_exclude = []
    # individual_only = False
    # path = '.'
    
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
    gps_dict = {}
    metadata_dict = {}
    
    errors = []
    warnings = []
    info = []
    
    #### Individual Metadata: 
    ### Initialize plots
    plt.close('all') #close any previous figures  
    
    ## Create a PDF output of plots
    plot_dir = 'plots'
    
            
    ##Create battery and temperature time series graphs for all SN
    batt_temp_fig = plt.figure(0)
    batt_temp_ax = batt_temp_fig.subplots(2)
    batt_temp_ax[0].set_title("Battery Voltage")
    batt_temp_ax[0].set_ylabel("voltage (V)")
    batt_temp_ax[0].set_xlabel("month-date hour")
    batt_temp_ax[1].set_title("Temperature")
    batt_temp_ax[1].set_ylabel("temperature (C)")
    batt_temp_ax[1].set_xlabel("month-date hour")
    batt_temp_fig.tight_layout()
    
    ##Create A2 and A3 time series graphs for all SN
    A2_A3_fig = plt.figure(1)
    A2_A3_ax = A2_A3_fig.subplots(2)
    A2_A3_ax[0].set_title("A2")
    A2_A3_ax[0].set_ylabel("Voltage (V)")
    A2_A3_ax[0].set_xlabel("month-date hour")
    A2_A3_ax[1].set_title("A3")
    A2_A3_ax[1].set_ylabel("Voltage (V)")
    A2_A3_ax[1].set_xlabel("month-date hour")
    A2_A3_fig.tight_layout()
    
    ##Create GPS runtime histogram plots for all SN 
    gps_fig = plt.figure(2)
    gps_ax = gps_fig.subplots(len(SN_list))
    gps_ax[0].set_title("GPS Runtime")
    gps_fig.tight_layout()
    
        ## Individual Metadata tests:
    for SN_index, SN in enumerate(SN_list):
        print('\nChecking metadata for ' + SN)
        metadata = pd.read_csv(path +'/metadata/' + SN + 'metadata_000.txt', sep = ',')
        
        #### battery voltage must be in reasonable range (1.7 to 15 V)
        pstats_df.loc[SN, "battery min"] = min(metadata.batt)
        pstats_df.loc[SN,"battery max"] = max(metadata.batt)
        #battery voltage minimum tests
        if pstats_df.loc[SN, "battery min"] < 1.7:
            errors_df.loc[SN, "battery min"] = "ERROR"
            err_message = f"{SN} BATTERY ERROR: battery level {np.round(1.7-min(metadata.batt),decimals=2)} Volts below minimum threshold (1.7V)."
            print(err_message)
            errors.append(err_message)
        elif pstats_df.loc[SN, "battery min"] < 3.0:
            errors_df.loc[SN, "battery min"] = "WARNING"
            warn_message = f"{SN} BATTERY WARNING: battery level within {np.round(min(metadata.batt-1.7),decimals=2)} Volts of minimum threshold (1.7V)."
            print(warn_message)
            warnings.append(warn_message)
        else:
            errors_df.loc[SN, "battery min"] = "OKAY"
            
        #battery voltage maximum tests    
        if pstats_df.loc[SN, "battery max"] > 15:
            errors_df.loc[SN, "battery max"] = "ERROR"
            err_message = f"{SN} BATTERY ERROR: battery level {np.round(max(metadata.batt)-15,decimals=2)} Volts above maximum threshold (15V)."
            print(err_message)
            errors.append(err_message)
        elif pstats_df.loc[SN, "battery max"] > 14.95:
            errors_df.loc[SN, "battery max"] = "WARNING"
            warn_message = f"{SN} BATTERY WARNING: battery level within {np.round(15 - max(metadata.batt-1.7), decimals=2)} Volts of maximum threshold (15V)."
            print(warn_message)
            warnings.append(warn_message)
        else:
            errors_df.loc[SN, "battery max"] = "OKAY"  
        
    ##%%%%%##
        ####temperature must be within reasonable range (-20 t0 60 C)    
        pstats_df.loc[SN, "temperature min"] = min(metadata.temp)
        pstats_df.loc[SN, "temperature max"] = max(metadata.temp)
        pstats_df.loc[SN, "temperature average"] = np.mean(metadata.temp)
     
        #temperature minimum check
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
        #temperature maximum check
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
        batt_temp_ax[0].plot(time_datestamp_dec, batt_dec)
        batt_temp_ax[0].legend(SN_list)
        
        #Plot temperature
        batt_temp_ax[1].plot(time_datestamp_dec, temp_dec)
        batt_temp_ax[1].legend(SN_list)
        
    
                  
    ##%%%%%##         
        #### A2 and A3 must be 0-3.1, and dV/dt = 0 should be true <1% of record
        #re-evaluate threshold %age
        #Add increased information in error messages
        
        ##A2
        A2_zerodiff_proportion = (np.sum(np.diff(metadata.A2) == 0) / (len(metadata.A2) -1 ))
        pstats_df.loc[SN, "A2 dV/dt nonzero"] = A2_zerodiff_proportion
        within_A2_range = (all(metadata.A2 >=0) & all(metadata.A2 <= 3.1))
        pstats_df.loc[SN, "A2 range"] = within_A2_range
        
        #Check that A2 dV/dt == 0 less than 99% of time. >99% indicates a likely short circuit.
        if pstats_df.loc[SN, "A2 dV/dt nonzero"] > 0.99: 
            errors_df.loc[SN, "A2 dV/dt nonzero"] = "ERROR"
            err_message = f"{SN} A2 ERROR: {np.round(A2_zerodiff_proportion*100,decimals=1)}% of A2 dV/dt is greater than 0. Less than 99% indicates a likely short circuit."
            print(err_message)
            errors.append(err_message)
        elif pstats_df.loc[SN,"A2 dV/dt nonzero"] > 0.95: #95% placeholder; uncertain interpretation
            errors_df.loc[SN, "A2 dV/dt nonzero"] = "WARNING"
            warn_message = f"{SN} A2 WARNING: {np.round(A2_zerodiff_proportion*100,decimals=1)}% of A2 dV/dt is greater than 0. Less than 99% indicates a likely short circuit."
            warnings.append(warn_message)
            print(warn_message)
        else:
            errors_df.loc[SN, "A2 dV/dt nonzero"] = "OKAY" 
      
        #failure_type = "A2"
        errors_df.loc[SN, "A2"] = np.NaN
        A2_check = (np.sum(np.diff(metadata.A2) == 0) / (len(metadata.A2) -1 ))
        limit = 0.01
        if A2_check > limit: 
            failure_message = SN + ': A2 dV/dt error'
            errors.append(failure_message)
            errors_df.loc[SN, "A2 dV/dt"] = A2_check
            print(f"{SN} A2 ERROR: A2 dV/dt constant ratio of {np.round(A2_check,decimals=3)} (limit is {limit}).")
        if not (all(metadata.A2 >=0) & all(metadata.A2 <= 3.1)):
            failure_message = SN + ': Bad A2'
            errors_df.loc[SN, "A2 error"] = 1 #yes error exists (true)
            errors.append(failure_message) 
            print(f"{SN} A2 ERROR: A2 malfunction")
        else:
            errors_df.loc[SN, "A2 error"] = 0 #no error present(false)
        errors_df.loc[SN,"A2 dV/dt nonzero"] = "OKAY"
            
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
        ##A3
        A3_nonzero = (np.sum(np.diff(metadata.A3) == 0) / (len(metadata.A3) -1 ))
        pstats_df.loc[SN,"A3 dV/dt nonzero"] = A3_nonzero
        within_A3_range = (all(metadata.A3 >=0) & all(metadata.A3 <= 3.1))
        pstats_df.loc[SN, "A3 range"] = within_A3_range
        
        #Check that A3 dV/dt == 0 less than 99% of time
        if pstats_df.loc[SN,"A3 dV/dt nonzero"] > 0.99: 
            errors_df.loc[SN,"A3 dV/dt nonzero"] = "ERROR"
            err_message = f"{SN} A3 ERROR: {np.round(A3_nonzero*100,decimals=1)}% of A3 dV/dt is greater than 0. Less than 99% indicates a likely short circuit."
            errors.append(err_message)
            print(err_message)
        elif pstats_df.loc[SN,"A3 dV/dt nonzero"] > 0.95: #placeholder of 95%
            errors_df.loc[SN, "A3 dV/dt nonzero"] = "WARNING"
            warn_message = f"{SN} A3 WARNING: {np.round(A3_nonzero*100,decimals=1)}% of A3 dV/dt is greater than 0. Less than 99% indicates a likely short circuit."
            warnings.append(warn_message)
            print(warn_message)
        else:
            errors_df.loc[SN,"A3 dV/dt nonzero"] = "OKAY"
            
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
        
        #Plot A3 data
        A2_A3_ax[1].plot(time_datestamp_dec, A3_dec)
        A2_A3_ax[1].legend(SN_list)
        
    ##%%%%%##            
        #### minFifoFree and maxFifoUsed should always add to 75
        fifo_sum = metadata.minFifoFree + metadata.maxFifoUsed
        pstats_df.loc[SN, "FIFO sum"] = max(fifo_sum)
        if any((fifo_sum) != 75):
           errors_df.loc[SN, "FIFO sum"] = "ERROR"
           err_message = f"{SN} FIFO ERROR: FIFO sum exceeds range by {np.round(75 - max(fifo_sum),decimals=2)}."
           errors.append(err_message)
           print(err_message)
        ##elif for warning??
        else:
            errors_df.loc[SN,"FIFO within range"] = "OKAY"
            
        #### maxFifoUsed should be less than 5 99% of the time, and should never exceed 25
        #how to display in dataframe if less than 5 99%?
        max_fifo_check = np.sum(metadata.maxFifoUsed > 5)/(len(metadata.maxFifoUsed) -1 )
        pstats_df.loc[SN,"max fifo within range"] = max_fifo_check
        if max_fifo_check > 0.01:
            errors_df.loc[SN,"max fifo within range"] = "INFO"
            info_message = f"{SN} MAX FIFO INFO: max fifo (max_fifo_check - 0.01)*100 % outside the standard operating range."
            print(info_message)
            info.append(info_message)
        else:
            errors_df.loc[SN, "max fifo within range"] = "OKAY"
            
        if any(metadata.maxFifoUsed > 25):
            errors_df.loc[SN, "max fifo"] = "WARNING"
            warn_message = f"{SN} MAX FIFO WARNING: max fifo exceeds maximum by {np.round(max(metadata.maxFifoUsed)-25,decimals=2)}."
            warnings.append(warn_message)
        else:
            errors_df.loc[SN, "max fifo within range"] = "OKAY"
            
    ##%%%%%##             
        #### maxOverruns should always be zero 
        pstats_df.loc[SN, "max overruns"] = max(metadata.maxOverruns)
        if any(metadata.maxOverruns) !=0:
            errors_df.loc[SN, "max overruns"] = "ERROR"
            err_message = f"{SN} OVERRUNS ERROR: maximum overruns does equal 0. ({max(metadata.maxOverruns)})"
            print(err_message)
            errors.append(err_message)
        else:
            errors_df.loc[SN, "max overruns"] = "OKAY"
            
    ##%%%%%##             
        #### unusedStack1 and unusedStackIdle should always be above some threshold 
        pstats_df.loc[SN,"unused stack1 max"] = max(metadata.unusedStack1)
        pstats_df.loc[SN,"unused stack idle max"] = max(metadata.unusedStackIdle)
        if any(metadata.unusedStack1 <= 30) or any(metadata.unusedStackIdle <= 30):
            errors_df.loc[SN, "unused stack"] = "WARNING"
            warn_message = f"UNUSED STACK WARNING: One value of unused stack exceeds maximum by {np.round(max(metadata.unusedStack1)-30,decimals=2)}."
            print(warn_message)
            warnings.append(warn_message)
        else:
            errors_df.loc[SN, "unused stack"] = "OKAY"

    ##%%%%%##             
        #### find time differences among samples with gps off that are > 180 sec
        time_check = np.diff(metadata.t[metadata.gpsOnFlag == 0])[10:] # skip the first ten seconds and first GPS cycle, which is very long by design
        pstats_df.loc[SN, "gps run time"] = max(time_check)
        if any(time_check > 180): 
            errors_df.loc[SN, "gps run time"] = "WARNING"
            warn_message = f"{SN} GPS WARNING: GPS runtime is {max(time_check)}."
            print(warn_message)
            warnings.append(warn_message)
        ## individual GPS:
            #plot GPS histogram for runtime
            #change plotting options 
        
        #IDEAS - change to stacked bar chart with different colors for serial numbers
        # how to select colors automatically?
        time_filt = time_check[time_check > 12]
        time_filt[time_filt > 180] = 180
        binsize = np.arange(10,180,10)
        gps_ax[SN_index].hist(time_filt, bins=np.arange(10,180,5))
        gps_ax[SN_index].set_ylabel('#' + SN_list[SN_index])
        gps_ax[SN_index].axes.xaxis.set_ticklabels([])
        gps_ax[SN_index].xaxis.set_major_locator(plt.MultipleLocator(20))
        gps_ax[SN_index].xaxis.set_minor_locator(plt.MultipleLocator(10))
        if SN == SN_list[-1]:    
            gps_ax[SN_index].set_xlabel('seconds')
            gps_ax[SN_index].axes.xaxis.set_ticklabels([0,20,40,60,80,100,120,140,160,'>180'])
        
           
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
    pd.set_option('display.max_colwidth', -1)
    print('Results')
    print(errors_df)
    print('Stats:')
    print(pstats_df)
    
    
    print("\nSerial number tests complete.") 
   
    s = pstats_df.shape
    print(s)
    ## Create a PDF output of plots
    report_date = datetime.datetime.today()
    report_date = report_date.strftime("%Y-%m-%d")
    filename = str("Huddle_test_output_" + report_date)
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font('helvetica', size=12)
    pdf.cell(0,10, f"Huddle Test Results: {report_date}", border=1, ln=0, align= 'C')
    pdf.ln() #new line
    pdf.cell(0,10, f"Date: {report_date}",border=1,align= 'C', ln=1)
    pdf.cell(50,10,"This is a column header", ln=1)
    pdf.set_font('helvetica', size=8)
    pdf_header = []
    #attempt to get multiline headers
    for index, col in enumerate(pstats_df.columns):
       col_header = [[],[]]
       col_header = col.split() 
       pdf_header.append(col_header)
    #replace temperature with temp
    for word in pdf_header:
        if word[0] == "temperature":
            word[0] = "temp"
    
    for i, header in enumerate(pdf_header):
        if len(pdf_header[i]) == 3:
            pdf_header[i].append('')
        elif len(pdf_header[i]) == 2:
            pdf_header[i].append('')
            pdf_header[i].append('')
        elif len(pdf_header[i]) == 1:
            pdf_header[i].append('')
            pdf_header[i].append('')
            pdf_header[i].append('')
        else:
            pass
    for i, header in enumerate(pdf_header):
        print(pdf_header)
        for j in range(0,4):
               pdf.set_font('helvetica', size=8)
               pdf.cell(12,5, '%s' % pdf_header[i][4])
       #pdf.ln() 
       
    pdf.ln()
    
    #print out stats data frame into pdf
    for i in range(0,len(pstats_df)):
       for j in range(0,len(pstats_df.columns)): 
           pdf.cell(12,10, '%s' % np.round((pstats_df.iloc[i,j]),3), 1, 0, 'C')
       pdf.ln()

    
    pdf.output(f"{filename}.pdf")
    #return {'errors':errors, 'warnings':warnings, 'stats':pstats_df, 'results':errors_df}
    
 #%%  
    ## Before running the group tests, ensure that we actually have data more than one Gem!
    ## If not, add a warning, and return without conducting group tests.
    if (len(SN_list) == 1) or individual_only:
        warnings.append('Test only includes one logger; cannot run comparison tests')
        return {'errors':errors, 'warnings':warnings, 'notes':notes}

    ## If we're at this point, we have data from multiple loggers and are clear to run group tests.
    ## Group metadata:
    #### all loggers' first and last times should agree within 20 minutes
    start_time = metadata_dict[SN_list[0]].t.min()
    stop_time = metadata_dict[SN_list[0]].t.max()
    failure_type = "time sync"
    for SN in SN_list[1:]:
        if np.abs(metadata_dict[SN].t.min() - start_time) > (20*60):
            failure_message = SN + ': metadata start times disagree excessively'
            print(failure_message)
            errors.append(failure_message)
        else:
            print(SN + ': metadata start times agree')
        if np.abs(metadata_dict[SN].t.max() - stop_time) > (20*60):
            failure_message = 'metadata stop times disagree excessively'
            print(failure_message)
            errors.append(failure_message)
        else:
            print(SN + ': metadata stop times agree')
        start_time = max(start_time, metadata_dict[SN].t.min())
        stop_time = min(stop_time, metadata_dict[SN].t.max())
        
    #### at every given time, temperature must agree within 2C for all loggers
    times_to_check = np.arange(start_time, stop_time)
    #average_temperatures = 0
    all_temperatures = np.zeros((len(SN_list), len(times_to_check)))
    temperatures = {}
    for i, SN in enumerate(SN_list):
        temperatures[SN] = scipy.signal.lfilter(np.ones(100)/100, [1],
            scipy.interpolate.interp1d(metadata_dict[SN].t, metadata_dict[SN].temp)(times_to_check))
        # to do: use the median instead
        #average_temperatures += temperatures[SN]/ len(SN_list)
        all_temperatures[i,:] =  temperatures[SN]
    med_temperatures = np.median(all_temperatures, 0)
    for SN in SN_list:
        if np.sum(np.abs(temperatures[SN] - med_temperatures) > 2)/len(times_to_check) > 0.1:
            failure_message = SN + ': disagrees excessively with average temperature'
            print(failure_message)
            errors.append(failure_message)
        else:
            print(SN + ': Temperatures agree')
    
    ## Group GPS
    #### all loggers' average lat and lon should agree within 1 m
    #### all loggers' first and last GPS times should agree within 20 minutes

    ## group waveform data data:
    #### length of converted data should match among all loggers
    #### a "coherent window" has all cross-correlation coefficients > 0.9, passes consistency criterion, and has amplitude above noise spec. 90% of coherent windows should have only nonzero lags, and none should have persistently nonzero lags (define).
    DB = gemlog.make_db(path + '/mseed', '*', 'tmp_db.csv')
    DB = DB.loc[DB.station.isin(SN_list),:]
    [t, lag, xc_coef, consistency] = check_lags(DB)
    coherent_windows = (consistency == 0) & (np.median(xc_coef, 0) > 0.8)
    zero_lags = lag[:,coherent_windows]==0
    num_coherent_windows = np.sum(coherent_windows)
    if (np.sum(np.all(zero_lags, 0)) / num_coherent_windows) < 0.8:
        failure_message = 'Time lags are excessively nonzero for coherent time windows'
        print(failure_message)
        errors.append(failure_message)
    else:
        print('Time lags for coherent time windows are mostly/all zero')
        

    return {'errors':errors, 'warnings':warnings, 'notes':notes}
