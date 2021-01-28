
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from scipy import signal
from scipy.signal import butter

#def butter_bandpass(lowcut, highcut, fs, order=5):
#    nyq = 0.5 * fs
#    low = lowcut / nyq
#    high = highcut / nyq
#    b, a = butter(order, [low, high], btype='band')
#    return b, a
#
#def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
#    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
#    y = lfilter(b, a, data)
#    return y
#
#def butter_lowpass(cutoff, fs, order=9):
#    nyq = 0.5 * fs
#    normal_cutoff = cutoff / nyq
#    b, a = butter(order, normal_cutoff, btype='low', analog=False)
#    return b, a
#
#def butter_lowpass_filter(data, cutoff, fs, order=9):
#    b, a = butter_lowpass(cutoff, fs, order=order)
#    y = lfilter(b, a, data)
#    return y
#
#def butter_highpass(cutoff, fs, order=3):
#    nyq = 0.5 * fs
#    normal_cutoff = cutoff / nyq
#    b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
#    return b, a
#
#def butter_highpass_filter(data, cutoff, fs, order=3):
#    b, a = butter_highpass(cutoff, fs, order=order)
#    y = signal.filtfilt(b, a, data)
#    return y

ds = xr.open_dataset("ao_weekly_1983-2018.nc")
#print(ds.time)
ds = ds.resample(time='QS-DEC').mean()
ds = ds.isel(time=slice(None,None,4))

print(ds)

z = ds.ao.values

print(z)

fig, ax = plt.subplots(1,1, figsize=(7,4))

ax.plot(ds.time, z, 'ko-')

ax.grid(True, color='grey', ls=':')

ax.set_title("Winter (DJF) mean AO index", fontsize=15)
ax.set_ylabel("Standardized value", fontsize=12)
ax.set_xlabel("Year", fontsize=12)

plt.show()

