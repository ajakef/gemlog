import numpy as np
from scipy.interpolate import CubicHermiteSpline, interp1d
from scipy.optimize import minimize
from scipy.stats import linregress
from scipy.ndimage import median_filter
import matplotlib.pyplot as plt
import gemlog
import pdb
#default_deg1 = 0.001024 # make this an argument; different for Gem vs Aspen 100/200 sps

def get_GPS_spline(G, default_deg1 = 0.001024):
    ## account for default slope
    x = np.array(G['msPPS'])# * default_deg1
    y = np.array(G['t'])# - x # in a perfect world, y is constant
    block_slope, block_mean_x, block_mean_y = get_block_results(x, y, default_deg1)

    ## comment out compatibility check. it's not wrong, but it's more aggressive than we need now.
    #for i in range(len(block_slope) - 1):
    #    if not check_compatibility(block_slope[i], block_mean_x[i], block_mean_y[i],
    #                               block_slope[i+i], block_mean_x[i+1], block_mean_y[i+1]):
    #        raise gemlog.core.CorruptRawFileDiscontinuousGPS
    
    #return CubicHermiteSpline(block_mean_x, block_mean_y, block_slope)
    return interp1d(block_mean_x, block_mean_y, fill_value = 'extrapolate')


## define blocks and estimate slopes and mean x/y values
def get_block_results(x, y, default_deg1):
    min_block_length = 4 # previously 8
    max_block_length = 30
    max_block_gap = 30/default_deg1
    max_block_sep = 1000/default_deg1 ## currently not used
    block_starts_all = [0]
    block_starts = []
    block_ends = []
    for i in range(1, len(x)):
        if ((x[i] - x[i-1]) > max_block_gap) or ((i - block_starts_all[-1]) > max_block_length):
            block_starts_all.append(i)
    block_starts_all.append(len(x))
    for i in range(0, len(block_starts_all) - 1):
        if (block_starts_all[i+1] - block_starts_all[i]) >= min_block_length:
            block_starts.append(block_starts_all[i])
            block_ends.append(block_starts_all[i+1])
    #block_starts.append(len(x))
    block_slope = np.zeros(len(block_starts))
    block_mean_x = np.zeros(len(block_starts))
    block_mean_y = np.zeros(len(block_starts))
    for i in range(len(block_starts)):
        xb = x[block_starts[i]:block_ends[i]]
        yb = y[block_starts[i]:block_ends[i]]
        block_slope[i], block_mean_x[i], block_mean_y[i] = _get_block_stats(xb, yb, default_deg1)
    if _check_step_between_blocks(block_mean_x, block_mean_y, block_slope, default_deg1)[0]:
        raise gemlog.exceptions.CorruptRawFileDiscontinuousGPS('Likely step between GPS cycles')
        
    # check for other timing errors between blocks (1/1800 is 1-sec step over 30 minutes)
    if any(np.abs(np.diff(block_mean_y)/np.diff(block_mean_x) - default_deg1)/default_deg1 > 1/1800):
        breakpoint()
        raise gemlog.exceptions.CorruptRawFileInadequateGPS('Timing error between GPS cycles')
    return block_slope, block_mean_x, block_mean_y

def _check_step_between_blocks(block_mean_x, block_mean_y, block_slope, default_deg1):
    # check for a step: good slopes, discordant intercepts
    step_x = []
    step_y = []
    for i in range(len(block_mean_x)-1):
        if np.abs(block_slope[i] - default_deg1)/default_deg1 < 0.01 and \
           np.abs(block_slope[i+1] - default_deg1)/default_deg1 < 0.01 and \
           np.abs(np.diff(block_mean_y[i:(i+2)]) - default_deg1 * np.diff(block_mean_x[i:(i+2)])) > 0.5:
            step_x.append(np.mean(block_mean_x[i:(i+2)]))
            step_y.append(np.mean(block_mean_y[i:(i+2)]))
    return (len(step_x) > 0), step_x, step_y

## check for slope-intercept compatibility between adjacent blocks
def check_compatibility(m1, x1, y1, m2, x2, y2, default_deg1):
    # find intersection x for lines defined by slope and point
    if (np.abs(y2-y1 + m2*(x1-x2)/2 - m1*(x2-x1)/2) < 0.0025) or \
       (np.abs(m1 - m2)/default_deg1 < 1e-5): # make sure halfway between them is returned in case they are collinear
        x = (x1 + x2)/2
    else:
        x = (y1 - y2 + m2*x2 - m1*x1) / (m2 - m1) # this works if they are not collinear (otherwise zero division can occur)

    # check that intersection is between the 2 blocks
    if (x < x1) or (x > x2):
        pdb.set_trace()
        return False
    # other checks? slopes close together? interval between x1 and x2 not too long?
    
    return True
    
    
    
##################################

# This estimator fits the best half of the data closely (L2) and totally ignores the other half.
# Note that this is NOT the same as L1-norm, which is affected by outliers but much less than L2.
# It also isn't trying to minimize the median absolute deviation, which is much less precise.
def _calc_abs_resid(a, b, x, y):
    return np.abs(y - a - b * (x - x[0]))

def _rms_sub_med(params, xy):
    x = xy[0]
    y = xy[1]
    resid = _calc_abs_resid(params[0], params[1], xy[0], xy[1])
    return np.sqrt(np.mean(resid[resid <= np.median(resid)]**2))

def _check_step_within_block(x, y, default_deg1):
    # this is a different problem from check_step_between_blocks because we don't have block_slope
    # and there might be glitches
    # detrend and use medfilt to ignore glitches and detect steps
    yfilt = median_filter(y - y[0] - (x-x[0]) * default_deg1, 3)
    w = np.where(np.abs(np.diff(yfilt)) > 0.5)[0]
    return len(w)>0, x[w], y[w]

    
def _get_block_stats(x, y, default_deg1, max_dev = 0.005, max_errors = 2):
    # Assume that GPS points are all either good, spikes, or steps.
    # If a brief (<len(x)/2) non-reversed step is present, raise an exception.
    # Reversed steps are treated the same as spikes. If either are present, drop them.
    if _check_step_within_block(x, y, default_deg1)[0]:
        raise gemlog.exceptions.CorruptRawFileDiscontinuousGPS('Excessive unfit points in GPS data, likely step')

    # at this point, we think there is not a step in this block. calculate the stats.
    # Calculate the line of best fit for the 50% of the data that can be fit best.
    # This completely ignores spikes and the short half of a step.
    # Default method BFGS has failed in tests ("success: False") so using Nelder-Mead.
    dydx_estimate = (y[-1] - y[0])/(x[-1] - x[0])
    result = minimize(_rms_sub_med, [y[0], dydx_estimate], [x, y], method = 'Nelder-Mead')
    if not result.success:
        #breakpoint()
        raise gemlog.core.CorruptRawFileInadequateGPS('Failed to fit GPS data')
    a, b = result.x # a intercept, b slope
        
    resid_excessive = _calc_abs_resid(a, b, x, y) > max_dev
    i_fit = np.where(~resid_excessive)[0]
    if np.sum(resid_excessive) > max_errors:
        # block couldn't be fit adequately as a line; check for a step to be safe
        result2 = minimize(_rms_sub_med, [y[0], dydx_estimate], [x[resid_excessive], y[resid_excessive]], method = 'Nelder-Mead')
        a2, b2 = result2.x
        resid_excessive2 = _calc_abs_resid(a2, b2, x[resid_excessive], y[resid_excessive]) > max_dev
        i_fit2 = np.where(resid_excessive)[0][~resid_excessive2]
        if result2.success and \
           (np.max(i_fit2) < np.min(i_fit)) or (np.max(i_fit) < np.min(i_fit2)) and \
           np.abs(b - default_deg1)/default_deg1 < 0.1 and \
           np.abs(b2 - default_deg1)/default_deg1 < 0.1 and \
           np.abs(a - a2) > 0.5:
            raise gemlog.exceptions.CorruptRawFileDiscontinuousGPS('Excessive unfit points in GPS data, likely step')
        else:
            raise gemlog.core.CorruptRawFileInadequateGPS('Excessive unfit points in GPS data')
    x = x[~resid_excessive]
    y = y[~resid_excessive]
    linreg = linregress(x, y)
    x_mean = (x[0] + x[-1])/2
    y_mean = (y[0] + y[-1])/2
    return linreg.slope, x_mean, y_mean
        
    
    
