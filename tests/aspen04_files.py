"""
Writes synthetic AspenCSV0.4 raw files the same way the Aspen firmware does, for tests. The formatting follows the
firmware's printHeader, printData, and printTimeAnchor (aspen-firmware/SD_functions.h) and formatTimeAnchorLine
(aspen-firmware/src/firmware_logic.cpp).
"""
import numpy as np


class Aspen04Writer:
    def __init__(self, filename, SN = '00022', gains = (128, 32, 32, 32), sample_int_ms = 5, use_compression = True):
        self.f = open(filename, 'w')
        self.gains = gains
        self.max_channel = max([i + 1 for i in range(4) if gains[i] != 0])
        self.use_compression = use_compression
        self.last = [0, 0, 0, 0]
        self.next_absolute = True
        self.f.write('#AspenCSV0.4\n')
        self.f.write('#ADC1,ADC2,ADC3,ADC4: data lines are exactly one sample interval (C line) apart; the line after each T line is absolute, the rest are deltas\n')
        self.f.write('#T,tcxoSec,utcSec,sync: time of the next data line (center of the ADC filter window); sync 0=no GPS time yet (utcSec blank), 1=GPS locked, 2=GPS holdover\n')
        self.f.write('#Trailing channels of a delta line with zero change (or gain 0) are omitted; first enabled channel is always present\n')
        self.f.write('#G,msPPS,msLag,yr,mo,day,hr,min,sec,lat,lon,z,HDOP,num_sats,age_ms\n')
        self.f.write('#H,ms,batt,V,mA,tempC,RH,Vcc,maxSampleDiff_ms,maxReadTime_us,maxWriteTime_us,maxDrainTime_us,drainDiff_ms,gpsOn,sync,timeErr_us,tcxo_ppm,adc_ppm,maxRetimeTime_us\n')
        self.f.write(f'S,{SN}\n')
        self.f.write('#Firmware0.8.0-0_2026-10-16T00:00:00Z\n')
        self.f.write('C,' + ','.join(str(g) for g in gains) + f',{sample_int_ms},1,15,25,0,0,1\n')

    def anchor(self, tcxo_s, utc_s, sync):
        """T line: the next data line is at this time, and is absolute."""
        tcxo_us = int(round(tcxo_s * 1e6))
        if sync == 0:
            self.f.write(f'T,{tcxo_us // 1000000}.{tcxo_us % 1000000:06d},,0\n')
        else:
            utc_ms = int(round(utc_s * 1000))
            self.f.write(f'T,{tcxo_us // 1000000}.{tcxo_us % 1000000:06d},{utc_ms // 1000}.{utc_ms % 1000:03d},{sync}\n')
        self.next_absolute = True

    def data(self, values):
        values = [int(v) for v in values] + [0] * (4 - len(values))
        absolute = self.next_absolute
        self.next_absolute = False
        deltas = [values[i] - self.last[i] for i in range(4)]
        compress = (not absolute) and self.use_compression and all(abs(d) <= 12 for d in deltas)
        if absolute:
            last_ch = self.max_channel - 1
        else:  # lastDeltaChannelIndex
            enabled = [i for i in range(4) if self.gains[i] != 0]
            last_ch = enabled[0]
            for i in enabled:
                if deltas[i] != 0:
                    last_ch = i
        chans = [i for i in range(last_ch + 1) if self.gains[i] != 0]
        if compress:
            self.f.write(''.join(chr(deltas[i] + 109) for i in chans) + '\n')
        else:
            self.f.write(','.join(str(values[i] if absolute else deltas[i]) for i in chans) + '\n')
        self.last = values

    def line(self, text):
        self.f.write(text + '\n')

    def close(self):
        self.f.close()


def health_line(ms, sync):
    return f'H,{ms},3.77,3.784,-8.0,22.7,31.9,3.29,6,1132,216,8375,101,1,{sync},-3.2,1.850,-12034,74'


def gps_line(utc_s, ms_pps):
    t = np.datetime64(int(utc_s), 's').astype(object)
    return (f'G,{ms_pps:.2f},120.00,{t.year},{t.month},{t.day},{t.hour},{t.minute},{t.second:.3f},'
            f'43.622260,-116.205940,820.0,0.95,9,350')
