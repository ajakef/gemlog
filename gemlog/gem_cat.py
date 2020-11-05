import getopt
import pdb
import warnings
import numpy as np
from numpy import NaN, Inf
import os, glob, csv, time, scipy
import pandas as pd
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
import obspy
import sys
import shutil
from gemlog.gemlog import _read_single, EmptyRawFile, CorruptRawFile

def gem_cat(input_dir, output_dir, ext = ''):
    """
    gem_cat
    Merge raw data files so that all contain GPS data.
    
    Translated from R code originally by Danny Bowman.
    
    Parameters
    ----------
    input_dir: raw gem directory
    output_dir: path for renumbered and concatenated files
    ext: extension of raw files to convert (normally the serial number; sometimes TXT for old Gems)
    
    Returns
    -------
    None; file output only
    """
    if not os.path.isdir(input_dir):
        raise Exception('Input path ' + output_dir + ' does not exist')
    
    if not os.path.isdir(output_dir):
        try:
            os.makedirs(output_dir)
        except:
            raise Exception('Output path ' + output_dir + ' does not exist and cannot be created')
        else:
            print('Created output folder ' + output_dir)
    gem_files = sorted(glob.glob(input_dir + '/' + 'FILE[0-9][0-9][0-9][0-9].' + ext + '*'))
    
    has_gps = np.zeros(len(gem_files))
    counter = 0
    out_file = [] ## note that out_file being empty acts as a flag saying that no files have been processed yet
    try:
        SN = pd.read_csv(gem_files[0], delimiter = ',', skiprows = 4, nrows=1, dtype = 'str', names=['s', 'SN']).SN[0]
    except:
        raise EmptyRawFile(gem_files[0])

    for k in range(len(gem_files)):
        #for k in range(2):        
        print(str(k+1) + ' of ' + str(len(gem_files)) + ': ' + gem_files[k])
        try:
            lines = pd.read_csv(gem_files[k], delimiter = '\n', dtype = 'str', names = ['line'])
        except:
            raise EmptyRawFile(gem_files[k])
        #gps_grep = grepl("^G.*-?\\d+\\.\\d+,-?\\d+\\.\\d+$", lines)
        gps_grep = lines.line.str.contains('^G')
        ## if no files have been processed yet and the current file has no GPS data, skip it
   
        if (sum(gps_grep) == 0) & (len(out_file) == 0):
            continue

        ## if this is the first file being processed, simply copy it and advance
        if len(out_file) == 0:
            #out_file = output_dir + "/FILE", sprintf("%04d", counter) + "." + SN
            out_file = output_dir + '/FILE%04d.%s' % (counter, SN)
            shutil.copy(gem_files[k], out_file)
            counter += 1
            has_gps[k] = 1
            continue

        if sum(gps_grep) > 0:
            has_gps[k] = 1
            if has_gps[k - 1] == 0:
                ## if this isn't the first file being processed and it does have GPS data but the previous file didn't, append it to the current outfile
                AppendFile(gem_files[k], out_file, gem_files[k-1])
            else:
                ## if this isn't the first file being processed and it has gps data and the previous file did too, start a new outfile
                #out_file = output_dir + "/FILE" + sprintf("%04d", counter) + "." + SN
                #file_copy(gem_files[k], out_file, overwrite = TRUE)
                out_file = output_dir + '/FILE%04d.%s' % (counter, SN)
                shutil.copy(gem_files[k], out_file)
                counter += 1
        else:
            ## if this isn't the first file being processed and there's no GPS data, append it to the current outfile
            #if k == 1:
            #    next
            #else:
                ##system(paste0("cat ", gem_files[k], " >> ",    out_file))
            AppendFile(gem_files[k], out_file, gem_files[k-1])
    return 

#%%
#infile = gem_files[1]
#outfile = out_file
#prev_infile = gem_files[0]
#if True:
def AppendFile(infile, outfile, prev_infile):
    # ensure that the output path exists
    outpath = os.path.dirname(outfile)
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        
    header = pd.read_csv(infile, delimiter = ',', nrows=1, dtype = 'str', names=['line']).line[0]
    format = header[7:]
    if format in ['0.85C', '0.9']:
        #SN = scan(infile, skip = 4, sep = ',', what = list(character(), character()), nlines = 1, flush = TRUE)[[1]][2]
        ## determine what the last ADC reading is before the start of this file
        #if len(prev_infile) == 12:
        #    path = '.'
        #else:
        #    path = prev_infile[0:-12]
        ###########################################
        #L = suppressWarnings(ReadGem(nums = num, path = path, units = 'counts')) # suppressWarnings because it'll warn that there's no GPS issue (which is kind of the point) or SN (which doesn't matter)
        #p_start = L$p[length(L$p)]
        p_start = int(_read_single(prev_infile, require_gps = False, version = format)['data'][-1,1])
        ########################################

        ## read the first data line of the infile and convert it to an ADC reading difference
        j=0
        while True:
            ###################
            linetype = pd.read_csv(infile, skiprows = j, delimiter = '\n', nrows=1, dtype = 'str', names=['s']).s[0][0]
            ###############
            if linetype == 'D':
                break
            j += 1

        ## adjust the data line and append it to the outfile
        s = pd.read_csv(infile, delimiter = ',', nrows=1, skiprows=j, dtype = 'str', names=['c1','c2'])
        with open(outfile, 'a') as OF, open(infile, 'r') as IF:
            OF.write(s['c1'].iloc[0] + ',' + str(int(s['c2'].iloc[0])-p_start) + '\n')
            ## append the rest of the file
            #system(paste0('tail -n +', j+2, ' ', infile, ' >> ', outfile))
            for k, line in enumerate(IF):
                if k >= (j+1):
                    OF.write(line)
                #if k > 10: # testing only
                #    break
    else:
        ## earlier format version
        ## Find the end of the header and append the infile to the outfile past that point
        j=0
        while True:
            linetype = pd.read_csv(infile, skiprows = j, delimiter = '\n', nrows=1, dtype = 'str', names=['s']).s[0][0]
            if linetype == 'D':
                break
            j += 1
        
        ## append the rest of the file
        #system(paste0('tail -n +', j+1, ' ', infile, ' >> ', outfile))
        with open(outfile, 'a') as OF, open(infile, 'r') as IF:
            #system(paste0('tail -n +', j+2, ' ', infile, ' >> ', outfile))
            for k, line in enumerate(IF):
                if k >= j:
                    OF.write(line)

    return

def print_call():
    print('gem_cat -i <input_dir> -o <output_dir> -e <ext>')
    print('input_dir: raw gem directory')
    print('output_dir: path for renumbered and concatenated files')
    print('ext: extension of raw files to convert (normally the serial number; sometimes TXT for old Gems)')
    

def main(argv = None):
    if argv is None:
        argv = sys.argv[1:]
    inputdir = 'raw'
    outputdir = 'raw_merged'
    ext = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:e:")
    except getopt.GetoptError:
        print_call()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print_call()
            sys.exit()
        elif opt in ("-i", "--inputdir"):
            inputdir = arg
        elif opt in ("-o", "--outputdir"):
            outputdir = arg
        elif opt in ("-e", "--ext"):
            ext = arg
    try:
        fn = os.listdir(inputdir)
    except:
        print("Problem opening input data folder " + inputdir + ". Did you give the right folder name after -i?")
        print("")
        print_call()
        sys.exit()
    if len(fn) == 0:
        print('Data folder ' + inputdir + ' does not contain any data files.')
        sys.exit()
    try:
        gem_cat(inputdir, outputdir, ext)
    except:
        print('gem_cat failed')
        sys.exit()
    else:
        print('Files merged successfully. Remember: sample times in the output are NOT precise!')
        

if __name__ == "__main__":
   main(sys.argv[1:])
