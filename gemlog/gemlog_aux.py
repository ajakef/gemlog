import matplotlib.pyplot as plt
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
