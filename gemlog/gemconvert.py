#/usr/bin/env python
# conda create -n g1 python=3.7.6 numpy obspy rpy2 pandas matplotlib
import sys # should always be available, doesn't need to be in "try"
try:
    import numpy as np
    import os, glob, getopt, logging, traceback, platform
    import gemlog
    from concurrent.futures import ProcessPoolExecutor
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

def convert_single_SN(arg_list):
    inputdir, SN, outputdir, output_format, output_length = arg_list
    logging.info(f'Beginning {SN}')
    try:
        #print([inputdir, SN, outputdir, output_format, output_length])
        gemlog.convert(inputdir, SN = SN, convertedpath = outputdir, output_format = output_format, file_length_hour = output_length)
        print(f'{SN} done')
    except KeyboardInterrupt:
        logging.info('Interrupted by user')
        print('Interrupted')
        #sys.exit()
    except Exception as e:
        print(e)
        ####logging.info(traceback.format_exc(e.__traceback__))
        logging.exception(parse_error(e), exc_info = gemlog._debug)
        ####logging.exception(traceback.format_exc())
        logging.info('Error in ' + SN)
        pass
    else:
        pass
    logging.info('Completed ' + SN)
    return 0

def print_call():
    print('gemconvert -i <inputdir> -s <serialnumbers> -x <exclude_serialnumbers> -o <outputdir> -f <format> -l <filelength_hours> -p <number_of_processes>')
    print('-i --inputdir: default ./raw/')
    print('-s --serialnumbers: separate by commas (no spaces); default all')
    print('-x --exclude_serialnumbers: separate by commas (no spaces); default none')
    print('-o --outputdir: default ./mseed')
    print('-f --format: mseed, sac, or tspair (text) currently supported; default mseed')
    print('-l --length: length of output converted files in hours; default 24')
    print('-t --test: if used, print the files to convert, but do not actually run conversion')
    print('-p --parallel: number of processes to run in parallel (limited by your computer); default 1.')
    print('-h --help: print this message')
    print('Problems: check/raise issues at https://github.com/ajakef/gemlog/issues/')
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
    output_length = 24 # hours
    num_processes = 1
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
        elif opt in ("-t", "--test"):
            test = True
        elif opt in ("-o", "--outputdir"):
            outputdir = arg
        elif opt in ("-f", "--format"):
            output_format = arg
        elif opt in ("-l", "--length"):
            output_length = arg
        elif opt in ("-p", "--parallel"):
            try:
                num_processes = int(arg)
            except:
                'could not interpret number_of_processes ' + arg + ' as integer'

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
    
    print(f'gemlog version {gemlog.__version__}')
    print('inputdir ', inputdir)
    print('serial numbers ', SN_list)
    print('outputdir ', outputdir)
    if not test:
        logging.info(f'***Starting conversion (gemlog version {gemlog.__version__})***')
        p = platform.uname()
        logging.info(f'Dependencies: Python {platform.python_version()}, obspy {gemlog.core.obspy.__version__}, pandas {gemlog.core.pd.__version__}, scipy {gemlog.core.scipy.__version__}, numpy {gemlog.core.np.__version__}')
        logging.info(f'platform.uname(): {p.system}, {p.release}, {p.version}, {p.machine}, {p.processor}')
        logging.info('Call: gemconvert ' + ' '.join(sys.argv[1:]))
        logging.info(f'inputdir="{inputdir}"')
        logging.info(f'outputdir="{outputdir}"')
        logging.info(f'serial number list = {SN_list}')
        logging.info(f'format="{output_format}", length_hours={output_length}, test={test}, parallel={num_processes}')
        
        #for SN in SN_list:
        with ProcessPoolExecutor(max_workers=num_processes) as pool:
            args_list = [[inputdir, SN, outputdir, output_format, output_length] for SN in SN_list]
            res = pool.map(convert_single_SN, args_list)

if __name__ == "__main__":
   main(sys.argv[1:])
