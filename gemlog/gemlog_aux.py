import matplotlib.pyplot as plt
import obspy
from gemlog import *
def PlotAmp(DB):
    allSta = DB.station.unique()
    allSta.sort()
    for sta in DB.station.unique():
        w = np.where(DB.station == sta)[0]
        w.sort()
        plt.plot(DB.t1[w], np.log10(DB.amp_HP[w]), '.')
        print(str(sta) + ' ' + str(np.quantile(DB.amp_HP[w], 0.25)))
    plt.legend(allSta)
    plt.show()

import obspy
def CheckDiscontinuity(files):
    st = obspy.Stream()
    for fn in files:
        st += obspy.read(fn)
    st.merge()
    st = st.split()
    if len(st) > 1:
        for tr in st:
            print(tr.stats.starttime.isoformat() + '--' + tr.stats.endtime.isoformat())


#import CheckDiscontinuity
import glob
fn = glob.glob('mseed/*095*')
CheckDiscontinuity(fn)


import sys
sys.path.append('/home/jake/Dropbox/Gem_logger/gem_package/python/gemlog_python/gemlog/gemlog')
import gemlog
gemlog.convert(rawpath = 'raw', SN = '095')
fn = glob.glob('raw/*')
fn.sort()
L=gemlog.ReadGem_v0_9(fn)
L2=gemlog.ReadGem(path='raw', SN ='095')


