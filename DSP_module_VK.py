# -*- coding: utf-8 -*-
"""
Created on Sat Aug  7 14:56:03 2021

@author: kundurv

"""

from scipy import signal
from scipy.fftpack import fft
import numpy as np
import statistics as stat

# %% Take FFT of the input data
def fft_blackmanharris(Din,time):
    """ computes FFT of input data with blackmanharris windowing.
        FFT length = input data length.
        Din = signal ( voltage, current etc).
        time = time values at which Din is recorded.
        time interval is expected to be constant.i.e t[1]-t[0] = t[n]-t[n-1] = Tsample.
    """
    N = len(Din)
    Tsample = time[1] - time[0];
    Fs = 1/Tsample;
    A_window = signal.blackmanharris(N)
    W1 = np.linalg.norm(A_window,1)
    W2 = np.linalg.norm(A_window,2)
    NRB = (W2/W1)**2
    NRBW =  Fs*NRB
    rms_w = np.sqrt((W2**2)/N)
    mean_w = stat.mean(A_window)
    ECF = 1/rms_w
    ACF = 1/mean_w
    DinT_window = np.multiply(A_window,Din)
    DoutFFT = np.abs(fft(DinT_window,N)*ECF*4/N).tolist() # FFT/N returns A/2 ( two sided spectrum). For single sided,*2 and for a fullscale sine FS, amplitude is FS/2, so *2 for dBFS
    DoutFFT_dB = 20*np.log10(DoutFFT)
    frequency = np.arange(Fs/N,Fs+Fs/N,Fs/N)
    return frequency, DoutFFT, DoutFFT_dB

# %% calculate SNDR of input FFT data
def SNDR_calc(DinF, Fsig1, Fsig2, Fs, BW_Fstart, BW_Fend):
    """ calculates the SNDR from the input fft data.
        DinF = FFT of signal. ( not in dB)
        Fsig1 = Signal frequency
        Fsig2 = Signal frequency. if single tone set this to 0.
        Fs = sampling frequency
        BW_Fstart = Bandwidth starting frequency
        BW_Fend = Bandwidth end frequency
    """
    N = int(len(DinF))
    Nstart = int(np.floor(N*BW_Fstart/Fs))
    Nend = int(np.ceil(N*BW_Fend/Fs))
    Nsig1 = int(np.ceil(N*Fsig1/Fs)) # signal bin
    Nsig2 = int(np.ceil(N*Fsig2/Fs)) # signal bin
    Nbin_sig = 6 # number of bins used for signal power calculation
    
    # %% signal power calculation
    Psig_watt = 0.0
    DinF_watt = (abs(np.array(DinF)))**2
    i = Nsig1 - np.floor(Nbin_sig/2)  # signal power is calculated over Nsig +/- Nbin_sig/2
    
    while i <=  Nsig1 + np.floor(Nbin_sig/2):
        Psig_watt = Psig_watt + DinF_watt[int(i)]
        i=i+1
    
    i =  Nsig2 - np.floor(Nbin_sig/2)
    if Fsig2 > 0 : # consider second tone only if Freq > 0
        while i <=  Nsig2 + np.floor(Nbin_sig/2):
            Psig_watt = Psig_watt + DinF_watt[int(i)]
            i=i+1
    
    #%% power at DC
    Pdc_watt = 0.0
    i=1
    while i <=  Nbin_sig:
        Pdc_watt = Pdc_watt + DinF_watt[int(i)]
        i=i+1
    
    #%% total power in BW
    Ptot_watt = 0.0
    i=Nstart
    while i <= Nend:
        Ptot_watt = Ptot_watt + DinF_watt[int(i)]
        i=i+1
    
    Pnoise_watt = Ptot_watt - Psig_watt
    
    #%% calculate noise in signal bins based on noise from average of neighboring  bins
    if Fsig2 > 0 :
        avg_noise_per_bin = Pnoise_watt/(Nend-Nstart-2*(Nbin_sig+1))
        Noise_in_sigbin = 2*(Nbin_sig+1)*avg_noise_per_bin
    else:
        avg_noise_per_bin = Pnoise_watt/(Nend-Nstart-2*(Nbin_sig+1))
        Noise_in_sigbin = (Nbin_sig+1)*avg_noise_per_bin
            
    Pnoise_watt = Pnoise_watt + Noise_in_sigbin
    
    SNDR = Psig_watt/Pnoise_watt
    SNDR_dB = 10*np.log10(SNDR)
    
    #%% Remove signal power, replace with avg noise power and compute spectrum
    
    Pfft_nosig = DinF[:]
    i = Nsig1 - np.floor(Nbin_sig/2)
    while i <= (Nsig1 + np.floor(Nbin_sig/2)):
        Pfft_nosig[int(i)] = np.sqrt(avg_noise_per_bin)
        i=i+1
    
    i = Nsig2 - np.floor(Nbin_sig/2)
    if Fsig2 > 0 :
        while i <= (Nsig2 + np.floor(Nbin_sig/2)):
            Pfft_nosig[int(i)] = np.sqrt(avg_noise_per_bin)
            i=i+1
    
    #%% calculate SFDR
    Pspur_max_inband = abs(max(Pfft_nosig[Nstart:Nend])) # Inband spur magnitude
    Pspur_max_OOBlo  = abs(max(Pfft_nosig[Nbin_sig:Nstart-1])) # Low side spur w.r.t desired band
    Pspur_max_OOBhi  = abs(max(Pfft_nosig[Nend+1:int(N/2)])) # High side spur w.r.t desired band
    Pspur_max_OOB = max(Pspur_max_OOBhi,Pspur_max_OOBlo)
    Psig_peak = max(DinF[Nsig1-int(Nbin_sig/2):Nsig1+int(Nbin_sig/2)])
    SFDR_dB_inband = 20*np.log10(Psig_peak/Pspur_max_inband)
    SFDR_dB_OOB = 20*np.log10(Psig_peak/Pspur_max_OOB)
    
    #%% calculate IM3 or HD3
    
    if Fsig2 > 0 :
        Nhd3_lo = int(np.ceil(N*(2*Fsig2-Fsig1)/Fs))
        Nhd3_hi = int(np.ceil(N*(2*Fsig1-Fsig2)/Fs))
        HD3 = max(Pfft_nosig[Nhd3_hi],Pfft_nosig[Nhd3_lo])
        HD3_dBc = 20*np.log10(HD3/Psig_peak)
        
    return  SNDR_dB, SFDR_dB_inband, SFDR_dB_OOB, HD3_dBc, Psig_watt, Pnoise_watt, Pfft_nosig

    
    
    
        