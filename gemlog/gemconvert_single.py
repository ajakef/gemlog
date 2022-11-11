#/usr/bin/env python
import sys # should always be available, doesn't need to be in "try"
try:
    import numpy as np
    import os, getopt, logging, platform
    import gemlog
    from concurrent.futures import ProcessPoolExecutor
except Exception as e:
    print('Either dependencies are missing, or the environment is not active')
    print('Error message:')
    print(e)
    sys.exit(2)

def print_call():
    print('gemconvert_single -i <input_file> -o <outputdir> -g <gps_required>')
    print('-i --input_file')
    print('-o --output_file: default ./<input_file>.mseed')
    print('-f --force-conversion: convert file even if it doesn\'t have enough GPS points for')
    print('accurate sample timing. Results will NOT be suitable for array processing')
    print('-h --help: print this message')
    print('Problems: check/raise issues at https://github.com/ajakef/gemlog/issues/')
    print('gemlog version: ' + gemlog.__version__)

def main(argv = None):
    if argv is None:
        argv = sys.argv[1:]
    input_file = ''
    output_file = None
    force_conversion = False
    try:
        opts, args = getopt.getopt(argv,"hfi:o:")
    except getopt.GetoptError:
        print_call()
        sys.exit(2)

    for opt, arg in opts:
        #print([opt, arg])
        if opt in ('-h', '--help'):
            print_call()
            sys.exit()
        elif opt in ("-i", "--input_file"):
            input_file = arg
        elif opt in ("-o", "--output_file"):
            output_file = arg
        elif opt in ("-f", "--force_conversion"):
            force_conversion = True

    print(f'gemlog version {gemlog.__version__}')
    print('input_file ', input_file)
    if force_conversion:
        print(' ')
        print('WARNING: Forcing conversion, even if the file doesn\'t have enough GPS data for')
        print('accurate sample timing. Do not use the resulting data for array processing or')
        print('anything else that requires accurate times.')
        print(' ')
    try:
        gemlog.core._convert_one_file(input_file, output_filename = output_file, require_gps = not force_conversion)
    except Exception as e:
        print(e)

if __name__ == "__main__":
   main(sys.argv[1:])
