import numpy as np; import gemlog;


## check that the cython reader works (most basic)
x = gemlog.gemlog._read_with_cython('FILE0110.110')

## 
d = gemlog.gemlog._read_single('./FILE0110.110')

gemlog.gemlog_aux._convert_raw_091_110('/home/jake/Dropbox/2022-01-23_SandiaSolarGems_4/raw/FILE0001.210', 'FILE0001.210')

## how fast does it run? faster than 0.91? similar? slower?

## does it give identical results? appears to so far.
