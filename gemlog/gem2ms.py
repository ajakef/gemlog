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
    print('gemconvert -i <inputdir> -s <serialnumbers> -x <exclude_serialnumbers> -o <outputdir> -f <format> -l <filelength_hours>')
    print('-i --inputdir: default ./raw/')
    print('-s --serialnumbers: separate by commas (no spaces); default all')
    print('-x --exclude_serialnumbers: separate by commas (no spaces); default none')
    print('-o --outputdir: default ./mseed')
    print('-f --format: mseed, sac, or tspair (text) currently supported; default mseed')
    print('-l --length: length of output converted files in hours; default 1')
    print('-t --test: if used, print the files to convert, but do not actually run conversion')
    print('-h --help: print this message')
    print('alias: gem2ms. gemlog version: ' + gemlog.__version__)

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
    #print(sys.executable)
    inputdir = 'raw'
    SN_list = ''
    exclude = []
    outputdir = None
    test = False
    output_format = 'MSEED'
    output_length = 1 # hours
    gemlog._debug = False
    try:
        opts, args = getopt.getopt(argv,"hdti:s:x:o:f:l:",["inputdir=","serialnumber="])
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
            output_format = arg
        elif opt in ("-l", "--length"):
            output_length = arg

    if outputdir is None:
        outputdir = output_format.lower()
    if gemlog._debug:
        logging.basicConfig(level=logging.DEBUG, filename="gemconvert_logfile.txt", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    else:
        logging.basicConfig(level=logging.INFO, filename="gemconvert_logfile.txt", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
        
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
        logging.info('***Starting conversion (gemlog version ' + gemlog.__version__ + ')***')
        logging.info('Call: gemconvert ' + ' '.join(sys.argv[1:]))
        for SN in SN_list:
            logging.info('Beginning ' + SN)
            try:
                gemlog.convert(inputdir, SN = SN, convertedpath = outputdir, output_format = output_format, file_length_hour = output_length)
            except KeyboardInterrupt:
                logging.info('Interrupted by user')
                print('Interrupted')
                sys.exit()
            except Exception as e:
                #logging.info(traceback.format_exc(e.__traceback__))
                logging.exception(parse_error(e), exc_info = gemlog._debug)
                #logging.exception(traceback.format_exc())
                logging.info('Error in ' + SN)
            else:
                logging.info('Completed ' + SN)

if __name__ == "__main__":
   main(sys.argv[1:])
