#/usr/bin/env python
import pandas as pd
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os, glob, getopt, sys, shutil
import obspy, gemlog
import time
import datetime
from gemlog.gemlog_aux import check_lags
from io import StringIO 
import pdb
from fpdf import FPDF
from matplotlib.backends.backend_pdf import PdfPages

def _metadata_status(status, message, status_list, serial_num, dataframe = None, col_name = None):
    # ordered input before named input
    """
    Writes the metadata status to a status dataframe and saves message to appropriate list.
    metadata_status(status, message, status_list, column_header, serial_num, dataframe = None)
    status : string - status of metadata (error, warning, note)
    message : string - resulting message that prints to console
    status_list : list - contains printed status messages
    col_name : string - name of metadata test
    serial_num : int - serial number of gem
    dataframe : pandas dataframe - Change from none if saving into individual test dataframe
    
    Returns
    -------
    status_list: updated status list with status message
    dataframe: and one word status to dataframe if specified

    """
    if dataframe is not None and col_name is not None:
        dataframe.loc[serial_num, col_name] = status.upper()
    status_list.append(message)
    print(message)
    return(dataframe, status_list)

def unique(list1):
    unique, index = np.unique(list1, return_index=True)
    return sorted(unique)

class PDF(FPDF):
    def footer(self):
        # Position at 1.5 cm from bottom
        self.set_y(-15)
        self.set_font('helvetica', 'I', 8)
        self.set_text_color(128)
        self.cell(0, 10, 'Page ' + str(self.page_no()), 0, 0, 'R')

    def heading(self, heading_name):
        # Arial 12
        self.ln()
        self.set_font('helvetica', 'B', size = 10)
        self.cell(0,10, heading_name, border=0, ln=0, align= 'L')
        self.ln()
        self.set_font('helvetica', size=8)
        
    def import_df(self, dataframe, SN_list):
        self.set_font('helvetica', size=6)
        for i in range(0,len(dataframe)):
           self.cell(8,10, '%s' % SN_list[i])
           for j in range(0,len(dataframe.columns)): 
               if type(dataframe.iloc[i,j]) is not str:
                   self.cell(12,10, '%s' % np.round(dataframe.iloc[i,j],3), 1, 0, 'C')
               else:
                   self.cell(12,10, '%s' % dataframe.iloc[i,j], 1, 0, 'C')
           self.ln()
           
    def import_list(self, list_name):
        self.set_font('helvetica', size=8)
        for i, value in enumerate(list_name):
            self.cell(12,4, '%s' % list_name[i], ln=1)
            
    def table_col(self, dataframe):
        status_header = []
        for index, col in enumerate(dataframe.columns):
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
            
        for j in range(0,4):
            if j == 2:
                self.cell(8,10, 'SN') # Create SN column heading
            else:
                self.cell(8,10, ' ') #Other lines are blank in SN column heading
            for i, header in enumerate(status_header):
                self.cell(12,5, '%s' % status_header[i][j])
            self.ln() 
            
def verify_huddle_test(path, SN_list = [], SN_to_exclude = [], individual_only = False, run_crosscorrelation_checks = False, generate_report = True):
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
    if False: ## set default input values in development; set to True if running the code line-by-line
        if os.getlogin() == 'tamara':
            #path = '/home/tamara/gemlog/demo_QC'
            path = '/home/tamara/gem_tests/huddle_test/2022-02-17_labtest_stop_time_bug'
        elif os.getlogin() == 'jake':
            path = '/home/jake/Work/gemlog_python/demo_QC'
        else:
            print('unknown user, need to define path')
        #SN_list = ['058','061','065','077']
        SN_list = []
        SN_to_exclude = []
        individual_only = False
        run_crosscorrelation_checks = False
        generate_report = True
        
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
    interval_dict = {}
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

    # Create wiggle figures
    wave_fig = plt.figure()
    wave_ax = wave_fig.subplots(len(SN_list))
    wave_ax[0].set_title('Waveforms')
    wave_fig.tight_layout()
    
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
        
        interval_dict[SN] = np.nanmean(np.diff(metadata.t)) # calculate interval between metadata sampling. use nanmean in case of possible missing times.
        
        ### Battery voltage must be in reasonable range
        # Define battery voltage range
        batt_min = 1.7
        batt_max = 15
        
        
        
        # Save a few statistics from battery metadata
        pstats_df.loc[SN, "battery min"] = min(metadata.batt)
        pstats_df.loc[SN,"battery max"] = max(metadata.batt)
       
        # Battery voltage minimum tests
        if pstats_df.loc[SN, "battery min"] < batt_min:
            err_message = f"{SN} BATTERY ERROR: battery level {np.round(1.7-min(metadata.batt),decimals=2)} Volts below minimum threshold (1.7V)."
            _metadata_status("error", err_message, errors, SN, dataframe=errors_df, col_name = "battery min")
        elif pstats_df.loc[SN, "battery min"] < batt_min + 1.5:
            warn_message = f"{SN} BATTERY WARNING: low battery level within {np.round(min(metadata.batt-1.7),decimals=2)} Volts of minimum threshold (1.7V)."
            _metadata_status("warning", warn_message, warnings, SN, dataframe=errors_df, col_name = "battery min")
        else:
            errors_df.loc[SN, "battery min"] = "OKAY"
            
        # Battery voltage maximum tests    
        if pstats_df.loc[SN, "battery max"] > batt_max:
            err_message = f"{SN} BATTERY ERROR: battery level {np.round(max(metadata.batt)-15,decimals=2)} Volts above maximum threshold (15V)."
            _metadata_status("error", err_message, errors, SN, dataframe=errors_df, col_name="battery max")
        elif pstats_df.loc[SN, "battery max"] > batt_max - 0.05:
            warn_message = f"{SN} BATTERY WARNING: battery level within {np.round(15 - max(metadata.batt-1.7), decimals=2)} Volts of maximum threshold (15V)."
            _metadata_status("warning", warn_message, warnings, SN, dataframe=errors_df, col_name="battery max")
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
            err_message = f"{SN} TEMPERATURE ERROR: temperature {np.abs(np.round(min(metadata.temp)+20,decimals=2))} degrees below minimum threshold (-20 C)."
            _metadata_status("error", err_message, errors, SN, dataframe=errors_df, col_name = "temperature min")
        elif pstats_df.loc[SN,"temperature min"] < -15: #modify as needed, just a backbone structure for now.
            warn_message = f"{SN} TEMPERATURE WARNING: temperature within {np.abs(np.round(20 + min(metadata.temp),decimals=2))} degrees of minimum threshold (-20 C)"
            _metadata_status("warning", warn_message, warnings, SN, dataframe=errors_df, col_name = "temperature min")
        else:
            errors_df.loc[SN, "temperature min"] = "OKAY" 
        # Temperature maximum check
        if pstats_df.loc[SN, "temperature max"] > 60: #degrees Celsius
            err_message = f"{SN} TEMPERATURE ERROR: temperature {np.round(max(metadata.temp)-60,decimals=2)} degrees above threshold (60 C)."
            _metadata_status("error", err_message, errors, SN, dataframe=errors_df, col_name = "temperature max")
        elif pstats_df.loc[SN, "temperature max"] > 50: #modify as needed, just a backbone structure for now.
            warn_message = f"{SN} TEMPERATURE WARNING: temperature within {np.round(60-max(metadata.temp),decimals=2)} degrees of maximum threshold (60 C)."
            _metadata_status("warning", warn_message, warnings, SN, dataframe=errors_df, col_name = "temperature max")
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
        # Horizontal line y-scales temperature too wide 
        #batt_temp_ax[1].axhline(temp_min, color="red", linestyle = ":", linewidth = 1)
        #batt_temp_ax[1].axhline(temp_max, color="red", linestyle=":", linewidth = 1)
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
            err_message = f"{SN} A2 ERROR: {np.round(A2_zerodiff_proportion*100,decimals=1)}% of A2 dV/dt is exactly 0. More than 99% indicates a likely short circuit."
            _metadata_status("error", err_message, errors, SN, dataframe=errors_df, col_name = "A2 flat")
        elif pstats_df.loc[SN,"A2 flat"] > 0.95: #95% placeholder; uncertain interpretation
            warn_message = f"{SN} A2 WARNING: {np.round(A2_zerodiff_proportion*100,decimals=1)}% of A2 dV/dt is exactly 0. More than 99% indicates a likely short circuit."
            _metadata_status("warning", warn_message, warnings, SN, dataframe=errors_df, col_name = "A2 flat")
        else:
            errors_df.loc[SN, "A2 flat"] = "OKAY" 
            
       #Check A2 is within range
        if not within_A2_range:
            pstats_df.loc[SN, "A2 range"] = 0 #outside of range (false)
            err_message = f"{SN} A2 ERROR: A2 outside of range"
            _metadata_status("error", err_message, errors, SN, dataframe=errors_df, col_name= "A2 range")
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
            err_message = f"{SN} A3 ERROR: {np.round(A3_nonzero*100,decimals=1)}% of A3 dV/dt is exactly 0. More than 99% indicates a likely short circuit."
            _metadata_status("error", err_message, errors, SN, dataframe=errors_df, col_name = "A3 flat")
        elif pstats_df.loc[SN,"A3 flat"] > 0.95: #placeholder of 95%
            warn_message = f"{SN} A3 WARNING: {np.round(A3_nonzero*100,decimals=1)}% of A3 dV/dt is exactly 0. More than 99% indicates a likely short circuit."
            _metadata_status("warning", warn_message, warnings, SN, dataframe=errors_df, col_name = "A3 flat")
        else:
            errors_df.loc[SN,"A3 flat"] = "OKAY"
            
       #Check A3 is within range     
        if not within_A3_range:
            pstats_df.loc[SN, "A3 range"] = 0 #outside of range (false)
            err_message = f"{SN} A3 ERROR: A3 outside of range."
            _metadata_status("error", err_message, errors, SN, dataframe=errors_df, col_name = "A3 range")
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
           notes_message = f"{SN} fifo notes: fifo sum exceeds range by {np.round(fifo_rule - max(fifo_sum),decimals=2)}."
           _metadata_status("note", notes_message, notes, SN, dataframe = errors_df, col_name = "fifo sum")
        ##elif for warning??
        else:
            errors_df.loc[SN,"fifo within range"] = "OKAY"
            
        #### maxFifoUsed should be less than 5 99% of the time, and should never exceed 25
        #how to display in dataframe if less than 5 99%?
        max_fifo_check = np.sum(metadata.maxFifoUsed > 5)/(len(metadata.maxFifoUsed) -1 )
        pstats_df.loc[SN,"max fifo within range"] = max_fifo_check
        if max_fifo_check > 0.01:
            notes_message = f"{SN} max fifo info: max fifo (max_fifo_check - 0.01)*100 % outside the standard operating range."
            _metadata_status("note", notes_message, notes, SN, dataframe = errors_df, col_name = "max fifo within range")
        else:
            errors_df.loc[SN, "max fifo within range"] = "OKAY"
            
        if any(metadata.maxFifoUsed > 25):  
            notes_message = f"{SN} max fifo note: max fifo exceeds maximum by {np.round(max(metadata.maxFifoUsed)-25,decimals=2)}."
            _metadata_status("note", notes_message, notes, SN, dataframe = errors_df, col_name = "max fifo")
        else:
            errors_df.loc[SN, "max fifo within range"] = "OKAY"
            
    ##%%%%%##             
        #### maxOverruns should always be zero 
        pstats_df.loc[SN, "max overruns"] = max(metadata.maxOverruns)
        if any(metadata.maxOverruns) !=0:
            warn_message = f"{SN} OVERRUNS WARNING: maximum overruns does not equal 0. ({max(metadata.maxOverruns)})"
            _metadata_status("warning", warn_message, warnings, SN, dataframe = errors_df, col_name = "max overruns")
        else:
            errors_df.loc[SN, "max overruns"] = "OKAY"
            
    ##%%%%%##             
        #### unusedStack1 and unusedStackIdle should always be above some threshold 
        #pstats_df.loc[SN,"unused stack1 max"] = max(metadata.unusedStack1)
        #pstats_df.loc[SN,"unused stack idle max"] = max(metadata.unusedStackIdle)
        if any(metadata.unusedStack1 <= 30) or any(metadata.unusedStackIdle <= 30):
            warn_message = f"UNUSED STACK WARNING: One value of unused stack exceeds maximum by {np.round(max(metadata.unusedStack1)-30,decimals=2)}."
            _metadata_status("warning", warn_message, warnings, SN, dataframe = errors_df, col_name = "unused stack")
        else:
            errors_df.loc[SN, "unused stack"] = "OKAY"

    ##%%%%%##             
        #### find time differences among samples with gps off that are > 180 sec
        # subtract interval from gps_time_check
        gps_time_check = np.diff(metadata.t[metadata.gpsOnFlag == 0])[10:] # skip the first ten seconds and first GPS cycle, which is very long by design
        gps_mean = gps_time_check[gps_time_check > 11] 
        pstats_df.loc[SN, "mean gps run time"] = np.mean(gps_mean)
        if any(gps_time_check > 180): 
            warn_message = f"{SN} GPS WARNING: GPS runtime is {max(gps_time_check)}, limit is 180."
            _metadata_status("warning", warn_message, warnings, SN, dataframe = errors_df, col_name = "gps run time")
            
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
        time_filt = gps_time_check[gps_time_check > 11] - interval_dict[SN]
        time_filt[time_filt > 180] = 180
        binsize = np.arange(10,180,5)
        bins = gps_ax[SN_index].hist(time_filt, bins=binsize)
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
        gps_fig.savefig(gps_fig_path, dpi=300, bbox_inches = 'tight', pad_inches = 0.1)
        
           
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
        # Define parameters
        stream = obspy.read(path +'/mseed/*..' + SN + '..HDF.mseed')
        stream.merge()
        # Check for clipping - if it is flatlined
        # filter then normalize right before plotting
        trace = stream[0]
        #trace.detrend()
        sps = trace.stats.sampling_rate
        npts = trace.stats.npts # save shortcut to number of points from stats inside trace
        wave_start = trace.stats.starttime
        wave_end = trace.stats.endtime

        
        d = 1/sps # seconds per sample (space between each sample in time)
        max_seconds = npts / sps # calcalate total number of seconds (maximum time value)
        t = np.arange(0, max_seconds, d) # create evenly spaced time values from 0 until the maximum to match to data
        
        trim_seconds = 60 * 5
        #trim_seconds = np.arange(0,trim_seconds, d)
        
        trace.trim(wave_start+trim_seconds, wave_end-trim_seconds)
        trace.detrend()
        trace.normalize()
        wave_ax[SN_index].plot(t[0:len(trace)],trace)
        wave_ax[SN_index].set_ylim(-1,1)
        wave_ax[SN_index].set_ylabel(SN)
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0.30)
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
    group_warn = []
    group_notes = []
    
    print("\n\nRunning group GPS tests")
    if (len(SN_list) == 1) or individual_only:
        warn_message = 'Test only includes one logger; cannot run comparison tests'
        _metadata_status("warning", warn_message, group_warn, SN)
        #return {'errors':errors, 'warnings':warnings, 'notes':notes}
        
#### all loggers' first and last times should agree within 20 minutes
    max_mins = 20
    max_secs = max_mins * 60
    # Record all start and stop times in group test dataframe
    for SN_index, SN in enumerate(SN_list):
        metadata = metadata_dict[SN]
        start_time = min(metadata.t)
        stop_time = max(metadata.t)
        group_df.loc[SN, "start time"] = start_time
        group_df.loc[SN, "end time"] = stop_time
    
    # Find the median and create a range within the max time distance for start and stop times
    
    gps_start_min = np.min(group_df.iloc[:,0]) # find the median of start times
    gps_stop_max = np.max(group_df.iloc[:,1]) # find the median of stop times
    upper_start = gps_start_min + max_secs # create an upper bound for start times
    lower_stop = gps_stop_max - max_secs # create a lower bound for stop times
    
    # Check all SN start and stop times to ensure they are within range
    for SN_index, SN in enumerate(SN_list):
        if not group_df.iloc[SN_index,0] <= upper_start: # create bound for lower start times
            err_message = (f"{SN} GROUP GPS ERROR: The start time is {np.round((group_df.iloc[SN_index,0]-upper_start)/60,2)} minutes after the acceptable start time range.")
            _metadata_status("error", err_message, group_err,SN)
        if not lower_stop <= group_df.iloc[SN_index,1]:
            err_message = (f"{SN} GROUP GPS ERROR: The stop time is {np.round((lower_stop - group_df.iloc[SN_index,1])/60,2)} minutes before the acceptable end time range.")
            _metadata_status("error", err_message, group_err, SN)
    if len(group_err) == 0:
        note = "The start and stop times for all loggers agree."
        _metadata_status("note", note, group_notes, SN)

#### at every given time, temperature must agree within 2C for all loggers
# interval of 061 and 065 is 10 seconds
    check_interval = 60 #seconds (time between checks)
    all_interval= interval_dict.values()
    df_width = check_interval / int(min(all_interval)) # use largest interval to calculate size
   
    
    # Determine start and stop times (Unix time)
    mod = gps_start_min % df_width # modulus remainder to round start time to even minute
    temp_start = gps_start_min - mod # start time on an even minute
    mod = gps_stop_max % df_width
    temp_end = gps_stop_max - mod - check_interval*10 # stop time on an even minute, ten minutes before end
    
    # Create dataframe to house temperatures to check
    column_index = np.arange(0,int((temp_end-temp_start)/df_width),1) # theoretically, the number of values between start and end
    
    group_temp_df = pd.DataFrame(index = SN_list, columns = column_index ) # create dataframe to house temperatures at each minute for each SN
    times_checked = np.zeros((1,max(column_index)))
    
    for SN_index, SN in enumerate(SN_list):
        metadata = metadata_dict[SN]
        temp_interval = check_interval / int(interval_dict[SN])
        argstart_list = [] # reset closest start list
        argend_list = [] # reset closest end list
        times = metadata.t # call all the time metadata
        for time in times: # for each time value in time metadata
            argstart_list.append(np.abs(temp_start - time)) # create a list of time values minus start time (find closest to start)
            argend_list.append(np.abs(temp_end - time)) # create list of time values minus end time (find closest index to end)
        
        start_index = np.argmin(argstart_list) # find index for closest minute time for start
        #format to not include nans from metadata
        stop_index = np.argmin(argend_list) # find index for last time value
        # too high

        times_to_check_index = np.arange(start_index, stop_index, temp_interval) # create evenly spaced array of even minutes
        times_checked = np.arange(temp_start, temp_end, temp_interval) # must create index based on mutally agreed start time
        for df_index, index in enumerate(times_to_check_index):
            #TROUBLESHOOT: 061 and 065 saving into dataframe as nan after index 29
            group_temp_df.iloc[SN_index,df_index] = metadata.temp[index] # save minute temperature reading into dataframe
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
            ts = int(times_checked[column])
            time_lookup = datetime.datetime.utcfromtimestamp(ts)
            err_message = (f"SN {x} recorded temperatures greater than 1 on either side of the temperature median {temp_median} on {time_lookup}. Total temperature range = {temp_range} (alpha)")
            _metadata_status("error", err_message, group_err, SN)
    if error == False:
        print("The recorded temperatures are within two degrees Celcius (alpha_version)")  

     
            
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
            err_message = 'Time lags are excessively nonzero for coherent time windows'
            # metadata_status() not used because this block does not use SN
            print(err_message)
            group_err.append(err_message)
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
    info.append(f"""BATTERY TEST INFO: This test is designed to ensure the voltage of each gem is within [{batt_min} to {batt_max} Volts at each recorded instance. Gems within this specified range, but within a range 0.5 Volts from the threshold will result in a warning message. This message is printed into the console and into the automatically generated pdf within the metadata folder. A plot is generated to visualize which serial numbers are malfunctioning and where 
they are operating. This plot is also saved in the metadata folder under an automatically generated folder called figures.""")
    trouble.append("""BATTERY TROUBLESHOOTING: If you received a battery warning, it is likely the gem was not able to record any waveform data. This can usually be fixed by changing the batteries. If you received a battery error [INSERT PROBLEM]
[INSERT TROUBLESHOOT]""")
                   
    ## Temperature test information and troubleshooting ##
    info.append(f"""TEMPERATURE TEST INFO: This test ensures the gemlogger is recording ambient air temperatures within an reasonable range. This range is set at {temp_min} to {temp_max} Celsius or {(temp_min * 9/5) + 32} to {(temp_max * 9/5) + 32} 
Fahrenheit. At temperatures outside this range, the electrical components of the gem could malfunction.""")
    trouble.append("""TEMPERATURE TROUBLESHOOTING: If you received a temperature warning, you are approaching the limit of the temperature range operation for the gemlogger ({temp_min} to {temp_max} C). If this value does not reflect an accurate 
ambient air temperature, [INSERT TROUBLESHOOTING]""")
                   
    ## A2 and A3 test information and troubleshooting ##
    info.append(f"""A2 AND A3 TEST INFO: These tests ensures that the A2 and A3 connections on the circuit board are functioning properly. Metadata values are recorded as Voltage (details). The first test ensures that the metadata has not flatlined 
for more than 99% of the recorded time. The second test ensures it is within a range of {A_min} - {A_max}. A plot is generated to visualize A2 and A3 voltages and is saved in the metadata folder under an automatically generated folder 
called figures.""")
    trouble.append(f""" A2 AND A3 TROUBLESHOOTING: Developers are still trying to create a practical range for these values. If you received an error for A2 or A3, most likely the sensor is fine. Cause for concern would be a flat line indicated on the 
                   A2 or A3 plots, or values that greatly exceed the limits. [INSERT MORE TROUBLESHOOTING]""")
    ## FIFO 
    info.append(f"""FIFO TEST INFO:  Under development""")
    trouble.append(f"""FIFO TROUBLESHOOTING:  Under development""") 

    ## Max Overruns
    info.append(f"""MAX OVERRUNS INFO: Under development""")     
    trouble.append(f"""MAX OVERRUNS TROUBLESHOOTING: Under development""")

    ##          
                 
#%%
   
# ============================================================================= 
    # Create a PDF output of plots with the date of report, errors warning and
    # notes list, and metadata summary dataframes
    # Create seperate package to reduce gemlog dependencies for detailed report
# =============================================================================
# To Do:   
#  - side by side images

    if generate_report == True:
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
        pdf = PDF()
        if len(pstats_df.columns) > 10 or len(errors_df.columns) > 10:
            pdf.add_page(orientation = 'L')
            p_width = 275
            # 25% larger images in landscape mode
            img_height = 150
            img_width = 220
            landscape = True
        else:
            pdf.add_page(orientation = 'P')
            p_width = 175
            img_height = 120
            img_width = 176
            landscape = False
        pdf.set_font('helvetica', 'B', size=12)
        pdf.cell(0,10, "Huddle Test Results", border=0, ln=1, align= 'C')
        pdf.cell(0,10, f"Date: {report_date}",border=0, ln=1, align= 'C')
    
    ## Individual test results
        # Insert error and warning list
        pdf.heading("Errors and Warnings")
        pdf.import_list(errors)
        pdf.import_list(warnings)
        pdf.ln()
        
        # Insert notes list       
        pdf.cell(0,5,"Notes", align = 'L', ln=1)
        pdf.import_list(notes)
        
        # Insert sensor status for each test from errors_df
        pdf.heading('Gem Sensor Status Summary')
        pdf.table_col(errors_df)
        pdf.import_df(errors_df, SN_list)
    
        # Insert statistics from pstats_df
        pdf.heading("Gem Sensor Metadata Details") # section heading
        pdf.table_col(pstats_df) # column headings
        pdf.import_df(pstats_df,SN_list) # table values   
        
        ## Add figures into report
        absc = pdf.get_x()
        ordin = pdf.get_y()
        pdf.ln()
        pdf.image(batt_temp_fig_path, w = img_width , h = img_height) 
        pdf.ln()
        pdf.image(A2_A3_fig_path, w = img_width, h = img_height)
        pdf.ln()
        pdf.image(gps_fig_path, w = img_width, h = 135)
        pdf.ln()
        
    ## Group test results
        pdf.heading("Group Test Results")
        pdf.import_list(group_err)    
        pdf.import_list(group_warn)
        pdf.import_list(group_notes)
        
    ## add the time lags if they were actually calculated
        if run_crosscorrelation_checks:
            pdf.ln()
            pdf.image(time_lags_fig_path, w = img_width, h = img_height)
            pdf.ln()    
    ## Test info and troubleshooting    
        pdf.heading('Test Info and Troubleshooting')
        for i, note in enumerate(info):
            pdf.multi_cell(p_width,5, '%s' %info[i])
            pdf.ln()
    ## Close and name file    
        pdf.output(f"{report_path}/{filename}.pdf")
        print("A pdf report has been generated")
    
##############################################
def print_call():
    print('gemconvert -i <inputdir> -s <serialnumbers> -x <exclude_serialnumbers>')
    print('-i --inputdir: default ./raw/')
    print('-s --serialnumbers: separate by commas (no spaces); default all')
    print('-x --exclude_serialnumbers: separate by commas (no spaces); default none')
    print('-h --help: print this message')
    print('gemlog version: ' + gemlog.__version__)

def main(argv = None):
    #verify_huddle_test(path, SN_list = [], SN_to_exclude = [], individual_only = False, run_crosscorrelation_checks = False, generate_report = True)
    if argv is None:
        argv = sys.argv[1:]
    #print(sys.executable)
    inputdir = '.'
    SN_list = ''
    exclude = []
    outputdir = None
    gemlog._debug = True
    try:
        opts, args = getopt.getopt(argv,"hdti:s:x:o:f:l:p:",["inputdir=","serialnumber="])
    except getopt.GetoptError:
        print_call()
        sys.exit(2)
    #print(2)
    #print(opts)
    for opt, arg in opts:
        #print([opt, arg])
        if opt in ('-h', '--help'):
            print_call()
            sys.exit()
        elif opt in ("-d", "--debug"):
            gemlog._debug = True
        elif opt in ("-i", "--inputdir"):
            inputdir = arg
        elif opt in ("-s", "--serialnumbers"):
            arg = arg.split(',')
            #print(arg)
            SN_list = arg
        elif opt in ("-x", "--exclude_serialnumbers"):
            arg = arg.split(',')
            #print(arg)
            exclude = arg
  
    verify_huddle_test(inputdir, SN_list, exclude)


## this goes last: executable terminal command
if __name__ == '__main__':
    main(sys.argv[1:])
