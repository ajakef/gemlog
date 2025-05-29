import numpy as np
from scipy.interpolate import CubicHermiteSpline
import matplotlib.pyplot as plt

def get_GPS_spline(G):
    ## account for default slope
    default_deg1 = 0.001024
    x = G['millis']# * default_deg1
    y = G['t'].t# - x # in a perfect world, y is constant
    block_slope, block_mean_x, block_mean_y = get_block_results(x, y)
    for i in range(len(block_slope) - 1):
        if not check_compatibility(block_slope[i], block_mean_x[i], block_mean_y[i],
                                   block_slope[i+i], block_mean_x[i+1], block_mean_y[i+1]):
            raise gemlog.core.CorruptRawFileDiscontinuousGPS
    
    return CubicHermiteSpline(block_mean_x, block_mean_y, block_slope)


## define blocks and estimate slopes and mean x/y values
def get_block_results(x, y):
    min_block_length = 10
    max_block_length = 30
    max_block_gap = 30
    max_block_sep = 1000
    block_starts_all = [0]
    block_starts = [0]
    for i in range(1, len(x)):
        if ((x[i] - x[i-1]) > max_block_gap) or ((i - block_starts[-1]) > max_block_length):
            block_starts_all.append(i)

    for i in range(1, len(block_starts_all)):
        if (block_starts_all[i] - block_starts_all[i-1]) >= min_block_length:
            block_starts.append(block_starts_all[i])

    block_slope = []
    block_mean_x = []
    block_mean_y = []
    for i in range(len(block_starts) - 1):
        xb = x[block_starts[i]:block_starts[i+1]]
        yb = y[block_starts[i]:block_starts[i+1]]
        block_slope[i], block_mean_x[i], block_mean_y[i] = _get_block_stats(xb, yb)
    return block_slope, block_mean_x, block_mean_y

## check for slope-intercept compatibility between adjacent blocks
def check_compatibility(m1, x1, y1, m2, x2, y2):
    # find intersection x for lines defined by slope and point
    if np.abs(y2-y1 + m2*(x1-x2)/2 - m1*(x2-x1)/2) < 0.0025 : # make sure halfway between them is returned in case they are collinear
        x = (x1 + x2)/2
    else:
        x = (x1 - y2 + m2*x2 - x1*x1) / (m2 - m1) # this works if they are not collinear

    # check that intersection is between the 2 blocks
    if (x < x1) or (x > x2):
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

def _get_block_stats(x, y, max_dev = 0.001, max_errors = 2):
    # Assume that GPS points are all either good, spikes, or steps.
    # If a brief (<len(x)/2) non-reversed step is present, raise an exception.
    # Reversed steps are treated the same as spikes. If either are present, drop them.

    # First, calculate the line of best fit for the 50% of the data that can be fit best.
    # This completely ignores spikes and the short half of a step.
    # Default method BFGS has failed in tests ("success: False") so using Nelder-Mead.
    dydx_estimate = (y[-1] - y[0])/(x[-1] - x[0])
    result = minimize(_rms_sub_med, [y[0], dydx_estimate], [x, y], method = 'Nelder-Mead')
    if not result.success:
        raise gemlog.core.CorruptRawFileInadequateGPS('Failed to fit GPS data')
    a, b = result.x
    resid_excessive = _calc_resid(a, b, x, y) > max_dev
    if np.sum(resid_excessive) > max_errors:
        raise gemlog.core.CorruptRawFileInadequateGPS('Excessive unfit points in GPS data, possible step')
    x = x[~resid_excessive]
    y = y[~resid_excessive]
    linreg = scipy.stats.linregress(x, y)
    x_mean = (x[0] + x[-1])/2
    y_mean = (y[0] + y[-1])/2
    return linreg.slope, x_mean, y_mean
        
    
    
