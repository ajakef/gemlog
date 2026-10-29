"""
Cython-based file parsers.
"""

import numpy as np
from libc.stdio cimport fopen, fclose, fgets, fseek, ftell, fgetc, rewind, FILE, sscanf, SEEK_END, SEEK_SET
from libc.string cimport strcpy, strcat

def parse_gemfile(filename, n_channels = 1, n_row = 1560000, dt_ms = 10, drop_bad_last_line = False, trailing_zeros_skipped = False):
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

    drop_bad_last_line : bool
        If True, ignore the last line of the file if it is incomplete (no terminating newline) or
        cannot be parsed. This can happen when power is lost while the logger is writing.
        Intended for AspenCSV0.2 and higher; the default (False) preserves the original behavior.

    trailing_zeros_skipped : bool
        If True, data lines omit trailing deltas of 0 (AspenCSV0.3).


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

    # detect an incomplete last line (no terminating newline) once, before parsing
    cdef bint file_incomplete = False
    cdef bint last_bad = False # last kept line failed to parse completely
    cdef bint last_skipped = False # last line read was not recognized, so not kept
    cdef bint drop_last = drop_bad_last_line
    cdef bint skip0 = trailing_zeros_skipped
    cdef int k = 0, start = 0, comma = 0
    cdef long file_size = 0
    if drop_last:
        fseek(cfile, 0, SEEK_END)
        file_size = ftell(cfile)
        if file_size > 0:
            fseek(cfile, file_size - 1, SEEK_SET)
            file_incomplete = (fgetc(cfile) != 10)
        rewind(cfile)

    cdef char line[256]
    cdef char* read
    cdef char line_type = 0

    cdef int n_matched = 0


    # D placeholders
    cdef double DmsSamp = 0
    cdef int offset = 0
    # build an array of args to sscanf in D lines (addresses to "values" array)
    cdef int ADC0 = 0
    cdef int ADC1 = 0
    cdef int ADC2 = 0
    cdef int ADC3 = 0
    cdef int j = 0
    cdef D_with_comma = False
    
    # Build the D line format string
    cdef char format[20]
    strcpy(format, b'D%lf')  # Start with the double (DmsSamp)
    for i in range(n_channels):
        strcat(format, ",%d")  # Append ",%d" for each integer column
    
    # G placeholders
    cdef int yr = 0, mo = 0, day = 0, hr = 0, mn = 0
    cdef double msPPS = 0, msLag = 0, sec = 0, lat = 0, lon = 0
    # M/H placeholders
    cdef int maxLag = 0, minFree = 0, maxUsed = 0, maxOver = 0
    cdef int gpsFlag = 0, freeStack1 = 0, freeStackIdle = 0
    cdef double ms = 0, batt = 0, temp = 0, A2 = 0, A3 = 0
    cdef double V = 0, mA = 0, RH = 0, Vcc = 0
    cdef int maxSampleDiff = 0, maxReadTime = 0, maxWriteTime = 0
    cdef int maxDrainTime = 0, drainDiff = 0, gpsOnFlag = 0

    # array to store parsed data
    result_array = np.zeros((n_row, 12), dtype=np.double)
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
    cdef double[:] millis_view = result_millis

    cdef Py_ssize_t line_number = 0
    while True:
        read = fgets(line, sizeof(line), cfile);
        if read == NULL:
            # EOF
            break
        #print(line)
        #print(n_channels)
        line_type = line[0]
        #print(line_type)

        ## AspenCSV0.3: trailing zero deltas are omitted. Only lines starting with D or z carry a time
        ## prefix (D<millis>[,deltas] or z[deltas], where z means dt+1); all others are pure deltas.
        if skip0 and ( ((line_type >= 97) and (line_type <= 122)) or line_type == 68 or line_type == 45 or
                       ((line_type >= 48) and (line_type <= 57)) ):
            for i in range(n_channels):
                view[line_number, i] = 0
            ADC0 = 0
            ADC1 = 0
            ADC2 = 0
            ADC3 = 0
            current_dD_millis = (prev_dD_millis + dt_ms) % (2**13)
            last_bad = False
            start = 0 # index where the deltas begin
            comma = 0 # nonzero if there is no delta field at all
            if line_type == 122:
                current_dD_millis = (prev_dD_millis + dt_ms + 1) % (2**13)
                start = 1
            elif line_type == 68:
                n_matched = sscanf(line + 1, "%lf", &DmsSamp)
                last_bad = n_matched < 1
                current_dD_millis = DmsSamp
                j = 1
                while (line[j] != 44) and (line[j] != 10) and (line[j] != 0):
                    j = j + 1
                if line[j] == 44:
                    start = j + 1
                else:
                    comma = 1
            if comma == 0:
                if (line[start] >= 97) and (line[start] <= 122): # letter-coded deltas
                    k = start
                    for i in range(n_channels):
                        if (line[k] < 97) or (line[k] > 122):
                            break
                        view[line_number, i] = line[k] - 109
                        k = k + 1
                elif (line[start] != 10) and (line[start] != 0): # comma-separated integers
                    n_matched = sscanf(line + start, "%d,%d,%d,%d", &ADC0, &ADC1, &ADC2, &ADC3)
                    if n_matched < 1:
                        last_bad = True
                    for i in range(n_channels):
                        if i < n_matched:
                            view[line_number, i] = [ADC0, ADC1, ADC2, ADC3][i]
            prev_dD_millis = current_dD_millis
            millis_view[line_number] = current_dD_millis
            line_type = 68

        ## Parse lines that are all-lowercase (fully compressed)
        elif ( ( (line[0] >= 97) and (line[0] <= 121) ) or  # e.g. q
           ( (line[0] == 122) and (line[1] >= 97) and (line[1] <= 122) ) ): # e.g. zm; ord('a', 'y', 'z')
            # if line is length n_channels, use default dt
            # example: q
            if not ((line[n_channels] >= 97) and (line[n_channels] <= 122)): # ord('a', 'z')
                offset = 0
                current_dD_millis = (prev_dD_millis + dt_ms) % (2**13)
            # if 1 more value than channels, first value is dt
            # example: zm
            else: 
                offset = 1
                if line[0] == 122: # ord('z'), code for dt being high by 1
                    current_dD_millis = (prev_dD_millis + dt_ms + 1) % (2**13) # diff_millis
                else:
                    current_dD_millis = (prev_dD_millis + dt_ms + line[0] - 109) % (2**13) # diff_millis
            for i in range(n_channels):
                view[line_number, i] = line[i+offset] - 109 # diff_ADC
            prev_dD_millis = current_dD_millis
            millis_view[line_number] = current_dD_millis
            last_bad = False
            line_type = 68 # ord('D') # because the D line is in an elif block, this is safe and the D line code won't be invoked

        ## Parse lines that are not all lowercase and begin with D,z,-,0-9 (partly compressed or not at all compressed)
        ## examples: D1012,-4; D4,h; -34; z28; 
        elif ( (line_type == 68) or (line_type == 122) or (line_type == 45) or # Dz- == 68,122,45
             ((line_type >= 48) and (line_type) <= 57) ):  # ord('0', '9') == 48,57
            #print(line_type)

            # Check if line starts with D and contains a comma
            D_with_comma = False
            if line_type == 68:
                j = 1
                while (line[j] != 10) and (line[j] != 0):  # until newline or null terminator
                    if line[j] == 44:  # comma (,)
                        D_with_comma = True
                        break
                    j = j + 1
            ## if the initial character is D or z, it won't be followed by a comma so skip it
            if line_type == 68:
                offset = 1
            elif line_type == 122:
                offset = 1
            else:
                offset = 0
            #print(offset)
            # try scanning it as if it's a line of comma-separated integers
            n_matched = sscanf(line + offset, "%lf,%d,%d,%d,%d", &DmsSamp, &ADC0, &ADC1, &ADC2, &ADC3)
            last_bad = n_matched < 1

            # if only one number was matched but there's a comma, it's of the form D1529,qmxk
            if (line_type == 68) and (n_matched == 1) and D_with_comma: 
                j = 1
                while (line[j] != 44) and (line[j] != 10):
                    j = j + 1
                for i in range(n_channels):
                    view[line_number, i] = line[j+i+1] - 109
            # otherwise, it's a list of comma-separated numbers, with or without millis diff
            elif n_matched == n_channels: # time count skipped in file
                for i in range(n_channels):
                    view[line_number, i] = [DmsSamp, ADC0, ADC1, ADC2][i]
                if line_type == 122:
                    DmsSamp = (prev_dD_millis + dt_ms + 1) % (2**13)
                else:
                    DmsSamp = (prev_dD_millis + dt_ms) % (2**13)
            else: # should be n_channels + 1 for the time count
                for i in range(n_channels):
                    #print((n_channels, line_number, i))
                    view[line_number, i] = [ADC0, ADC1, ADC2, ADC3][i]
            millis_view[line_number] = DmsSamp
            prev_dD_millis = DmsSamp
            line_type = 68

        ## Done with data lines; moving on to GPS lines
        elif line_type == 71:  # ord('G') == 71
            # G,msPPS,msLag,yr,mo,day,hr,min,sec,lat,lon
            # G,8171,70,2020,6,20,5,21,22.0,43.62226,-116.20594
            n_matched = sscanf(line + 2,
                               "%lf,%lf,%d,%d,%d,%d,%d,%lf,%lf,%lf",
                               &msPPS, &msLag, &yr, &mo, &day, &hr, &mn,
                               &sec, &lat, &lon)
            last_bad = n_matched != 10
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

        elif line_type == 72:  # ord('H') == 72
            # H,millis,batt,V,mA,temp,RH,Vcc,maxSampleDiff,maxReadTime,maxWriteTime,maxDrainTime,drainDiff,gpsOnFlag
            # H,773,3.77,3.784,-8.0,22.7,31.9,3.29,6,1132,216,8375,101,0
            n_matched = sscanf(line + 2,
                               "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d,%d,%d,%d,%d",
                               &ms, &batt, &V, &mA, &temp, &RH, &Vcc,
                               &maxSampleDiff, &maxReadTime, &maxWriteTime,
                               &maxDrainTime, &drainDiff, &gpsOnFlag)
            last_bad = n_matched != 13
            view[line_number, 0] = batt
            view[line_number, 1] = V
            view[line_number, 2] = mA
            view[line_number, 3] = temp
            view[line_number, 4] = RH
            view[line_number, 5] = Vcc
            view[line_number, 6] = maxSampleDiff
            view[line_number, 7] = maxReadTime
            view[line_number, 8] = maxWriteTime
            view[line_number, 9] = maxDrainTime
            view[line_number, 10] = drainDiff
            view[line_number, 11] = gpsOnFlag
            millis_view[line_number] = ms

        elif line_type == 77:  # ord('M') == 77
            # M,ms,batt(V),temp(C),A2,A3,maxLag,minFree,maxUsed,maxOver,
            # gpsFlag,freeStack1,freeStackIdle
            # M,8001,3.02,22.1,1.412,2.052,94,66,9,0,0,57,86
            n_matched = sscanf(line + 2,
                               "%lf,%lf,%lf,%lf,%lf,%d,%d,%d,%d,%d,%d,%d",
                               &ms, &batt, &temp, &A2, &A3, &maxLag,
                               &minFree, &maxUsed, &maxOver, &gpsFlag,
                               &freeStack1, &freeStackIdle)
            last_bad = n_matched != 12
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
            last_skipped = True
            continue

        last_skipped = False
        type_view[line_number] = line_type
        line_number += 1

    fclose(cfile)

    # an unrecognized trailing line was never kept, so there is nothing to drop
    if drop_last and line_number > 0 and not last_skipped and (file_incomplete or last_bad):
        line_number -= 1

    return (
        result_array[:line_number, :],
        result_linetypes[:line_number],
        result_millis[:line_number],
    )


def parse_aspen04(filename, n_channels = 1, n_row = 1560000, drop_bad_last_line = True):
    """
    Cythonized parser for AspenCSV0.4 files.

    In this format, data lines carry no time. They are exactly one sample interval apart, and each T line gives
    the time of the data line after it. That data line holds absolute ADC values; all other data lines hold deltas
    from the previous data line. Trailing zero deltas may be omitted (they are 0).

    Parameters
    ----------
    filename : bytes
        The filename to parse. Must be of type `bytes` -- use filename.encode('utf-8') if needed.

    n_channels : int
        Number of channels being logged.

    n_row : int
        Max number of data lines to be able to store (the number of lines in the file is always enough).

    drop_bad_last_line : bool
        If True, ignore the last line of the file if it is incomplete (no terminating newline) or a data line that
        cannot be parsed. This can happen when power is lost while the logger is writing.

    Returns
    -------
    tuple:
        - 2-d array (number of data lines, n_channels) of the values exactly as written (absolute or delta)
        - list of (data_index, line) for each T, G, and H line, where data_index is the number of data lines before
          it and line is the raw line (bytes)
    """
    cdef char* fname = filename

    cdef FILE* cfile
    cfile = fopen(fname, "rb")
    if cfile == NULL:
        msg = "No such file or directory: '{}'".format(filename)
        raise FileNotFoundError(2, msg)

    # detect an incomplete last line (no terminating newline) once, before parsing
    cdef bint file_incomplete = False
    cdef long file_size = 0
    if drop_bad_last_line:
        fseek(cfile, 0, SEEK_END)
        file_size = ftell(cfile)
        if file_size > 0:
            fseek(cfile, file_size - 1, SEEK_SET)
            file_incomplete = (fgetc(cfile) != 10)
        rewind(cfile)

    cdef char line[256]
    cdef char* read
    cdef char c = 0
    cdef int n_matched = 0, i = 0, k = 0
    cdef int adc[4]
    cdef int last_kind = 0 # kind of the last line read: 0 skipped, 1 data, 2 T/G/H
    cdef bint last_bad = False # last data line failed to parse completely

    result_array = np.zeros((n_row, n_channels), dtype=np.double)
    cdef double[:, :] view = result_array
    others = []
    cdef Py_ssize_t n = 0

    while True:
        read = fgets(line, sizeof(line), cfile)
        if read == NULL:
            break # EOF
        c = line[0]
        ## data lines: letter-coded deltas (a-y) or comma-separated integers (start with a digit or '-')
        if ((c >= 97) and (c <= 121)) or (c == 45) or ((c >= 48) and (c <= 57)):
            if n >= n_row:
                fclose(cfile)
                raise ValueError('parse_aspen04: more data lines than n_row')
            last_bad = False
            if c >= 97:
                k = 0
                for i in range(n_channels):
                    if (line[k] < 97) or (line[k] > 121):
                        break
                    view[n, i] = line[k] - 109
                    k = k + 1
                if (line[k] != 10) and (line[k] != 13) and (line[k] != 0):
                    last_bad = True
            else:
                n_matched = sscanf(line, "%d,%d,%d,%d", &adc[0], &adc[1], &adc[2], &adc[3])
                if n_matched < 1:
                    last_bad = True
                for i in range(n_channels):
                    if i < n_matched:
                        view[n, i] = adc[i]
            n += 1
            last_kind = 1
        ## T, G, and H lines are rare: keep the raw text and parse it in python
        elif (c == 84) or (c == 71) or (c == 72): # ord('T'), ord('G'), ord('H')
            others.append((n, <bytes>line))
            last_kind = 2
        else:
            last_kind = 0

    fclose(cfile)

    if drop_bad_last_line:
        if (last_kind == 1) and (file_incomplete or last_bad):
            n -= 1
        elif (last_kind == 2) and file_incomplete:
            others.pop()

    return result_array[:n, :], others
