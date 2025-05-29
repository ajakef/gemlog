import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize
from scipy.stats import linregress
import matplotlib.pyplot as plt
import gemlog
import pdb
default_deg1 = 0.001024

def get_GPS_spline(G):
    ## account for default slope
    x = np.array(G['msPPS'])# * default_deg1
    y = np.array(G['t'])# - x # in a perfect world, y is constant
    block_slope, block_mean_x, block_mean_y = get_block_results(x, y)

    ## comment out compatibility check. it's not wrong, but it's more aggressive than we need now.
    #for i in range(len(block_slope) - 1):
    #    if not check_compatibility(block_slope[i], block_mean_x[i], block_mean_y[i],
    #                               block_slope[i+i], block_mean_x[i+1], block_mean_y[i+1]):
    #        raise gemlog.core.CorruptRawFileDiscontinuousGPS
    
    return CubicHermiteSpline(block_mean_x, block_mean_y, block_slope)


## define blocks and estimate slopes and mean x/y values
def get_block_results(x, y):
    min_block_length = 8
    max_block_length = 30
    max_block_gap = 30/default_deg1
    max_block_sep = 1000/default_deg1
    block_starts_all = [0]
    block_starts = []
    for i in range(1, len(x)):
        if ((x[i] - x[i-1]) > max_block_gap) or ((i - block_starts_all[-1]) > max_block_length):
            block_starts_all.append(i)
    block_starts_all.append(len(x))
    for i in range(0, len(block_starts_all) - 1):
        if (block_starts_all[i+1] - block_starts_all[i]) >= min_block_length:
            block_starts.append(block_starts_all[i])
    block_starts.append(len(x))
    block_slope = np.zeros(len(block_starts) - 1)
    block_mean_x = np.zeros(len(block_starts) - 1)
    block_mean_y = np.zeros(len(block_starts) - 1)
    for i in range(len(block_starts) - 1):
        xb = x[block_starts[i]:block_starts[i+1]]
        yb = y[block_starts[i]:block_starts[i+1]]
        block_slope[i], block_mean_x[i], block_mean_y[i] = _get_block_stats(xb, yb)
    return block_slope, block_mean_x, block_mean_y

## check for slope-intercept compatibility between adjacent blocks
def check_compatibility(m1, x1, y1, m2, x2, y2):
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
def _calc_resid(a, b, x, y):
    return y - a - b * (x - x[0])

def _rms_sub_med(params, xy):
    x = xy[0]
    y = xy[1]
    resid = _calc_resid(params[0], params[1], xy[0], xy[1])
    return np.sqrt(np.mean(resid[resid <= np.median(resid)]**2))

def _get_block_stats(x, y, max_dev = 0.005, max_errors = 2):
    # Assume that GPS points are all either good, spikes, or steps.
    # If a brief (<len(x)/2) non-reversed step is present, raise an exception.
    # Reversed steps are treated the same as spikes. If either are present, drop them.

    # First, calculate the line of best fit for the 50% of the data that can be fit best.
    # This completely ignores spikes and the short half of a step.
    # Default method BFGS has failed in tests ("success: False") so using Nelder-Mead.
    dydx_estimate = (y[-1] - y[0])/(x[-1] - x[0])
    result = minimize(_rms_sub_med, [y[0], dydx_estimate], [x, y], method = 'Nelder-Mead')
    if not result.success:
        breakpoint()
        raise gemlog.core.CorruptRawFileInadequateGPS('Failed to fit GPS data')
    a, b = result.x
    resid_excessive = _calc_resid(a, b, x, y) > max_dev
    if np.sum(resid_excessive) > max_errors:
        breakpoint()
        raise gemlog.core.CorruptRawFileInadequateGPS('Excessive unfit points in GPS data, possible step')
    x = x[~resid_excessive]
    y = y[~resid_excessive]
    linreg = linregress(x, y)
    x_mean = (x[0] + x[-1])/2
    y_mean = (y[0] + y[-1])/2
    return linreg.slope, x_mean, y_mean
        
    
    
