"""
Cython-based file parsers.
"""

import numpy as np
from libc.stdio cimport fopen, fclose, fgets, FILE, sscanf
from libc.string cimport strcpy, strcat

def parse_gemfile(filename, n_channels = 1, n_row = 1560000, dt_ms = 10):
    """
    Cythonized gem logfile parser.

    Parameters
    ----------
    filename : bytes
        The filename to parse. Must be of type `bytes` -- use
        filename.encode('utf-8') if needed.
	
    n_channels : int
        Number of channels being logged (1 for Gem, up to 4 for Aspen).
	
    n_row : int
        Max number of rows to be able to store. 1560000 is enough for all normal
        Gem files. Set higher when concatenating many files.

    dt_ms : int
        Nominal sample interval in milliseconds (10 for Gem, 5 for Aspen at 200 sps)


    Returns
    -------
    tuple of three aligned numpy arrays:

        - 2-d array of numeric values read from the file
        - 1-d array of characters indicating the line type
        - 1-d array of the millisecond value of the row
    """
    cdef char* fname = filename

    cdef FILE* cfile
    cfile = fopen(fname, "rb")
    if cfile == NULL:
        msg = "No such file or directory: '{}'".format(filename)
        raise FileNotFoundError(2, msg)

    cdef char line[256]
    cdef char* read
    cdef char line_type = 0

    cdef int n_matched = 0


    # D placeholders
    #cdef int ADC = 0
    #ADC_values = np.zeros(n_channels, dtype = int)
    cdef int ADC_values[4]
    cdef double DmsSamp = 0

    # Build the D line format string
    cdef char format[20]
    strcpy(format, b'D%lf')  # Start with the double (DmsSamp)
    for i in range(n_channels):
        strcat(format, ",%d")  # Append ",%d" for each integer column

    # build an array of args to sscanf in D lines (addresses to "values" array)
    cdef int ADC0 = 0
    cdef int ADC1 = 0
    cdef int ADC2 = 0
    cdef int ADC3 = 0
    
    # G placeholders
    cdef int yr = 0, mo = 0, day = 0, hr = 0, mn = 0
    cdef double msPPS = 0, msLag = 0, sec = 0, lat = 0, lon = 0
    # M placeholders
    cdef int maxLag = 0, minFree = 0, maxUsed = 0, maxOver = 0
    cdef int gpsFlag = 0, freeStack1 = 0, freeStackIdle = 0
    cdef double ms = 0, batt = 0, temp = 0, A2 = 0, A3 = 0

    # array to store parsed data
    #n_row = 780000  # max number of rows to expect: 750000 + 15000 + 15000
    #n_row = 1560000  # max number of rows to expect: 2*(750000 + 15000 + 15000)
    result_array = np.zeros((n_row, 11), dtype=np.double)
    # make a view for faster indexing.
    # see https://cython.readthedocs.io/en/latest/src/userguide/numpy_tutorial.html#efficient-indexing-with-memoryviews
    cdef double[:, :] view = result_array
    cdef double prev_dD_millis = 0, current_dD_millis = 0

    # 1-D array to store linetype (single chars)
    result_linetypes = np.zeros(n_row, dtype='c')
    cdef char[:] type_view = result_linetypes
    # 1-D array to store millis.
    # range is 0 to 2**13, so choose short int
    result_millis = np.zeros(n_row, dtype=np.double)
    # cdef short[:] millis_view = result_millis
    cdef double[:] millis_view = result_millis

    cdef Py_ssize_t line_number = 0
    # were this python 3.8 we could maybe use the walrus operator.  alas
    while True:
        read = fgets(line, sizeof(line), cfile);
        if read == NULL:
            # EOF
            break

        line_type = line[0]
        if (line_type >= 97) and (line_type <= 122): # ord('a'), ord('z')
            if (line[1] < 97) or(line[1] > 122):
                view[line_number, 0] = line[0] - 109 # diff_ADC
                current_dD_millis = (prev_dD_millis + dt_ms) % (2**13)
            else:
                current_dD_millis = (prev_dD_millis + dt_ms + line[0] - 109) % (2**13) # diff_millis
                view[line_number, 0] = line[1] - 109 # diff_ADC
            prev_dD_millis = current_dD_millis
            millis_view[line_number] = current_dD_millis
            line_type = 68 # ord('D') # because the D line is in an elif block, this is safe and the D line code won't be invoked
	    
        elif line_type == 68:  # ord('D') == 68
            # DmsSamp,ADC
            # D7780,-1
            #n_matched = sscanf(line + 1, "%lf,%d", &DmsSamp, &ADC)
            # Prepare the argument list
            # build an array of args to sscanf in D lines (addresses to "values" array)
            # Call sscanf with the dynamic format string and arguments
            n_matched = sscanf(line + 1, '%lf,%d,%d,%d,%d', &DmsSamp, &ADC0, &ADC1, &ADC2, &ADC3)

            if n_matched == n_channels: # time count skipped in file
                view[line_number, 0] = DmsSamp 
                view[line_number, 1] = ADC0
                view[line_number, 2] = ADC1
                view[line_number, 3] = ADC2
                DmsSamp = (prev_dD_millis + dt_ms) % (2**13)
            else: # should be n_channels + 1 for the time count
                view[line_number, 0] = ADC0
                view[line_number, 1] = ADC1
                view[line_number, 2] = ADC2
                view[line_number, 3] = ADC3
            millis_view[line_number] = DmsSamp
            prev_dD_millis = DmsSamp
                

        elif line_type == 71:  # ord('G') == 71
            # G,msPPS,msLag,yr,mo,day,hr,min,sec,lat,lon
            # G,8171,70,2020,6,20,5,21,22.0,43.62226,-116.20594
            n_matched = sscanf(line + 2,
                               "%lf,%lf,%d,%d,%d,%d,%d,%lf,%lf,%lf",
                               &msPPS, &msLag, &yr, &mo, &day, &hr, &mn,
                               &sec, &lat, &lon)
            view[line_number, 0] = msLag
            view[line_number, 1] = yr
            view[line_number, 2] = mo
            view[line_number, 3] = day
            view[line_number, 4] = hr
            view[line_number, 5] = mn
            view[line_number, 6] = sec
            view[line_number, 7] = lat
            view[line_number, 8] = lon
            millis_view[line_number] = msPPS

        elif (line_type == 77) or (line_type == 72):  # ord('M') == 77, ord('H') == 72
            # M,ms,batt(V),temp(C),A2,A3,maxLag,minFree,maxUsed,maxOver,
            # gpsFlag,freeStack1,freeStackIdle
            # M,8001,3.02,22.1,1.412,2.052,94,66,9,0,0,57,86
            n_matched = sscanf(line + 2,
                               "%lf,%lf,%lf,%lf,%lf,%d,%d,%d,%d,%d,%d,%d",
                               &ms, &batt, &temp, &A2, &A3, &maxLag,
                               &minFree, &maxUsed, &maxOver, &gpsFlag,
                               &freeStack1, &freeStackIdle)
            view[line_number, 0] = batt
            view[line_number, 1] = temp
            view[line_number, 2] = A2
            view[line_number, 3] = A3
            view[line_number, 4] = maxLag
            view[line_number, 5] = minFree
            view[line_number, 6] = maxUsed
            view[line_number, 7] = maxOver
            view[line_number, 8] = gpsFlag
            view[line_number, 9] = freeStack1
            view[line_number, 10] = freeStackIdle
            millis_view[line_number] = ms
        else:
            continue

        type_view[line_number] = line_type
        line_number += 1

    fclose(cfile)

    return (
        result_array[:line_number, :],
        result_linetypes[:line_number],
        result_millis[:line_number],
    )
