from gemlog.gemlog_aux import gem_noise, ims_noise
from gemlog.gem_network import get_gem_response, deconvolve_gem_response
import numpy as np
import pytest, shutil, os, obspy

def test_gem_noise_integral():
    noise_info = gem_noise(version = '1.0', freq_min = 0.1, freq_max = 20, spectype = 'amp')
    assert np.round(noise_info['noise'], 4) == 0.0039 # spec mentioned in datasheet
    
def test_IMS_noise_integral():
    noise_info = ims_noise('min', freq_min = 0.17, freq_max = 0.19, spectype = 'power')
    
    assert np.round(noise_info['noise'] - (0.0004 * 0.02), 6) == 0 ## spectral density of microbarom peak
    

def test_get_gem_response():
    get_gem_response(gain = 'high')
    get_gem_response(gain = 'low')

def test_deconvolve_gem_response():
    tr = obspy.read()[0]
    tr_corrected = deconvolve_gem_response(tr)
    assert np.round(tr_corrected.std() - 0.9989, 3) == 0
