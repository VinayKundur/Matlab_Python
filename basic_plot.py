import csv
import matplotlib.pyplot as plt
import numpy as np
from DSP_module_VK import fft_blackmanharris, SNDR_calc



# %% inputs
filepath = "C:\\Users\\kundurv\\Documents\\Projects-AIS\\Drone\\Drone_RX_40nm\\MATLAB\\CSV_files\\Async\\"
filename = "ADConly_m95dBm_Pblk48dBm10MHz_VCMBUF_Jph25ps_SS_Tsweep.csv"
fullscale = 1023 # Full scale of the signal. 10  bit ADC data
col = 2
Fsig1 = 4.4123e6 # frequency of first signal
Fsig2 = 3.5103e6 # frequency of second signal
BW_Fstart = 3.41e6 # bandwidth start frequency
BW_Fend = 4.51e6 # bandwidth end frequency


# %% Read csv
file = filepath + filename
common_mode = int(fullscale/2)
time = []
voltage = []

with open (file) as csv_file:
	csv_reader = csv.reader(csv_file, delimiter=',')
	line_count = 0
	for row in csv_reader:
		if line_count == 0:
			pass
		else:
			time.append(float(row[0]))
			voltage.append(float(row[col]))
		line_count = line_count + 1


Fsample = 1/(time[1]-time[0])
# %% plot input waveforms
plt.rcParams.update({'font.size': 16})

plt.figure(1)
plt.subplot(2,1,1)
plt.plot(np.array(time)*(1e6),voltage, color ='k',linewidth=2,label='ADC code')
plt.xlabel('time(us)')
plt.ylabel('ADC code')
# plt.title('ADC output')
plt.legend(loc = 'upper right')
plt.grid(True)

plt.subplot(2,1,2)
Nstart = int(len(time)/2) # plot 10% of points in the middle of the array
Npoints = 100
plt.plot(np.array(time[Nstart:Nstart+Npoints])*(1e6),voltage[Nstart:Nstart+Npoints], color ='k',linewidth=2,label='ADC code zoomed')
plt.xlabel('time(us)')
plt.ylabel('ADC code')
# plt.title('ADC output')
plt.legend(loc = 'upper right')
plt.grid(True)


# %% compute FFT and plot PSD
frequency,DoutFFT,DoutFFT_dB = fft_blackmanharris(np.array(voltage)-common_mode,time)
DoutFFT_dBFS = DoutFFT_dB - 20*np.log10(fullscale)
RBW = frequency[3]-frequency[2]
SNDR_dB, SFDR_dB_inband, SFDR_dB_OOB, HD3_dBc, Psig_watt, Pnoise_watt, DoutFFTnosig = SNDR_calc(DoutFFT, Fsig1, Fsig2, Fsample, BW_Fstart, BW_Fend)
DoutFFT_dB_nosig_dBFS = 20*np.log10(np.array(DoutFFTnosig)) - 20*np.log10(fullscale)

plt.figure(2)
plt.plot(frequency/1e6,DoutFFT_dBFS,color = 'k',linewidth=2, label='with signal')
plt.plot(frequency/1e6,DoutFFT_dB_nosig_dBFS,color = 'r',linewidth=2, label='without signal')
plt.ylabel('PSD(dBFS)')
plt.ylim([-100, 0])
plt.yticks(np.arange(-100,10,10))
plt.xlabel('Frequency(MHz)(RBW=' + str(round(RBW/1e3,3)) + 'KHz)')
plt.xlim([0,(Fsample/2)/1e6])
plt.xticks(np.arange(0,(Fsample/2)/1e6 + 1,1))
plt.grid(True)
plt.legend(loc = 'upper right')
plt.show()