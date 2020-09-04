#/usr/bin/env python
# conda create -n g1 python=3.7.6 numpy obspy rpy2 pandas matplotlib
import sys # should always be available, doesn't need to be in "try"
try:
    import numpy as np
    import os, glob, getopt, requests
    import logging, traceback
    import gemlog
except Exception as e:
    print('Either dependencies are missing, or the environment is not active')
    print('Error message:')
    print(e)
    sys.exit(2)

#%%
def find_SN(x):
    return x[9:13]

# function to get unique values from list while preserving order (set-based shortcut doesn't do this)
def old_unique(list1):  # O(n^2), bad
    # intilize a null list 
    unique_list = [] 
    # traverse for all elements 
    for x in list1: 
        # check if exists in unique_list or not 
        if x not in unique_list: 
            unique_list.append(x) 
    return unique_list

def unique(list1): # thanks Kevin!
    unique, index = np.unique(list1, return_index=True)
    return sorted(unique)

def print_call():
    print('gem2ms.py -i <inputdir> -s <serialnumbers> -x <exclude_serialnumbers> -o <outputdir> -f <format>')
    print('-i --inputdir: default ./raw/')
    print('-s --serialnumbers: separate by commas (no spaces); default all')
    print('-x --exclude_serialnumbers: separate by commas (no spaces); default none')
    print('-o --outputdir: default ./mseed')
    print('-f --format: mseed, sac, or wav; default mseed')
    print('-t --test: if used, print the files to convert, but do not actually run conversion')
    print('-h --help: print this message')
    print('gemlog version: ' + gemlog.__version__)

def parse_error(e):
    e = str(e)
    if(e.find('Unable to allocate') > 0):
        return 'Tried to allocate very large data array, probably due to large time gap between adjacent files'
    elif(e.find('NULLType') > 0):
        return 'Suspected malformed file'
    else:
        return e
    

def main(argv = None):
    if argv is None:
        argv = sys.argv[1:]
    if gemlog._debug:
        logging.basicConfig(level=logging.DEBUG, filename="gem2ms_logfile.txt", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    else:
        logging.basicConfig(level=logging.INFO, filename="gem2ms_logfile.txt", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
        
    #print(sys.executable)
    inputdir = 'raw'
    SN_list = ''
    exclude = []
    outputdir = 'mseed'
    test = False
    fmt = 'MSEED'
    gemlog._debug = False
    try:
        opts, args = getopt.getopt(argv,"hdti:s:x:o:f:",["inputdir=","serialnumber="])
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
        elif opt in ("-t", "--test"):
            test = True
        elif opt in ("-o", "--outputdir"):
            outputdir = arg
        elif opt in ("-f", "--format"):
            fmt = arg
            
    try:
        fn = os.listdir(inputdir)
    except:
        print("Problem opening input data folder " + inputdir + ". Did you give the right folder name after -i?")
        print("")
        print_call()
        sys.exit()

    if(len(SN_list) == 0): # if user does not provide SN_list, take unique SNs in order
        SN_list = unique([find_SN(fn[i]) for i in range(len(fn))]) # set takes unique values
        SN_list.sort()
    else: # if user provided SNs, keep the order, but take unique values
        SN_list = unique(SN_list)

    if (len(SN_list) == 0) & (len(opts) == 0):
        print("No data to convert; printing help message instead")
        print("")
        print_call()
        sys.exit()

    SN_list = [i for i in SN_list if i not in exclude]
    
    print('inputdir ', inputdir)
    print('serial numbers ', SN_list)
    print('outputdir ', outputdir)
    if not test:
        for SN in SN_list:
            logging.info('Beginning ' + SN)
            try:
                gemlog.Convert(inputdir, SN = SN, convertedpath = outputdir, fmt = fmt)
            except Exception as e:
                #logging.info(traceback.format_exc(e.__traceback__))
                logging.exception(parse_error(e))
                #logging.exception(traceback.format_exc())
                logging.info('Error in ' + SN)
            else:
                logging.info('Completed ' + SN)

if __name__ == "__main__":
   main(sys.argv[1:])

## errors:
## empty files on 043: rpy2.rinterface_lib.embedded.RRuntimeError: Error in if (t2 <= t1) { : missing value where TRUE/FALSE needed

## ./gem2ms.py -i ~/Work/Chile2020/2020-01-19_Summit/raw
## note that 018 mostly lacked GPS signal
## same error with 067 in ~/Work/Chile2020/2020-01-19_RadialNetwork
#[1] "File 12 of 12 : /home/jake/Work/Chile2020/2020-01-19_Summit/raw/FILE0035.018"
#[1] 720000 720000
#Traceback (most recent call last):
#  File "./gem2ms.py", line 54, in <module>
#    main(sys.argv[1:])
#  File "./gem2ms.py", line 51, in main
#    gemlog.Convert(inputdir, SN = SN, convertedpath = outputdir)
#  File "/home/jake/Dropbox/Gem_logger/gem_package/python/gemlog.py", line 169, in Convert
#    L = ReadGemPy(nums[(nums >= n1) & (nums < (n1 + (12*blockdays)))], rawpath, alloutput = False, requireGPS = True, SN = SN, units = 'counts', network = network, station = station, location = location)
#  File "/home/jake/Dropbox/Gem_logger/gem_package/python/gemlog.py", line 254, in ReadGemPy
#    gps = pd.DataFrame.from_dict({ key : np.asarray(L[4].rx2(key)) for key in L[4].names })
#TypeError: 'NULLType' object is not iterable


#~/Work/Chile2020/2020-01-19_RadialNetwork, 017
#[1] "File 11 of 11 : raw//FILE0097.017"
#[1] 720000 720000
#R[write to console]: Error in seq.int(0, to0 - from, by) : wrong sign in 'by' argument
#R[write to console]: In addition: 
#R[write to console]: There were 50 or more warnings (use warnings() to see the first 50)
#R[write to console]: 
#Traceback (most recent call last):
#  File "/home/jake/Dropbox/Gem_logger/gem_package/python/gem2ms.py", line 78, in <module>
#    main(sys.argv[1:])
#  File "/home/jake/Dropbox/Gem_logger/gem_package/python/gem2ms.py", line 75, in main
#     gemlog.Convert(inputdir, SN = SN, convertedpath = outputdir)
#   File "/home/jake/Dropbox/Gem_logger/gem_package/python/gemlog.py", line 169, in Convert
#     L = ReadGemPy(nums[(nums >= n1) & (nums < (n1 + (12*blockdays)))], rawpath, alloutput = False, requireGPS = True, SN = SN, units = 'counts', network = network, station = station, location = location)
#   File "/home/jake/Dropbox/Gem_logger/gem_package/python/gemlog.py", line 252, in ReadGemPy
#     LI = gemlog.InterpTime(L)
#   File "/home/jake/.local/lib/python3.6/site-packages/rpy2/robjects/functions.py", line 192, in __call__
#     .__call__(*args, **kwargs))
#   File "/home/jake/.local/lib/python3.6/site-packages/rpy2/robjects/functions.py", line 121, in __call__
#     res = super(Function, self).__call__(*new_args, **new_kwargs)
#   File "/home/jake/.local/lib/python3.6/site-packages/rpy2/rinterface_lib/conversion.py", line 40, in _
#     cdata = function(*args, **kwargs)
#   File "/home/jake/.local/lib/python3.6/site-packages/rpy2/rinterface.py", line 789, in __call__
#     raise embedded.RRuntimeError(_rinterface._geterrmessage())
# rpy2.rinterface_lib.embedded.RRuntimeError: Error in seq.int(0, to0 - from, by) : wrong sign in 'by' argument


