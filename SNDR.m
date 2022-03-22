%% Author: kundurv <kundurv@SCLKUNDURVLT01>
% Created: 2020-03-30
% SNDR calculation

%% INPUTS
%Dinf = FFT data
%Fsig1 = Signal Frequency first tone
%Fsig2 = Signal Frequency second tone. set this 0 if there is only one tone
%signal
%Fs=sampling frequency
%Fstart = starting frequency of the BW
%Fend = ending frequency of the BW

%% OUTPUTS
%SNDR_dB = SNDR calculated over (Fstart-Fend)
% P_sig_watt = signal power
% P_noise_watt = noise power

%% Script

function [SNDR_dB, SFDR_dB_inband,SFDR_dB_OOB,HD3_dBc, P_sig_watt, P_noise_watt , Pfft_nosig] = SNDR (Dinf, Fsig1,Fsig2, Fs, Fstart,Fend)
N=length(Dinf);
Nstart = ceil(N*Fstart/Fs)-1; % BW starting frequency
Nend = ceil(N*Fend/Fs)+1; % BW ending frequency
Nsig1 = ceil(N*Fsig1/Fs); % Signal Bin
Nsig2 = ceil(N*Fsig2/Fs); % Signal Bin
Nbin_sig = 6 ; % number of bins used for Signal power calculation

%% calculate signal power
P_sig_watt = 0;

i= Nsig1 - floor(Nbin_sig/2); % signal power is calculated over Nbin +/- Nbin_sign/2
while i <= (Nsig1 + floor(Nbin_sig/2))  
    P_sig_watt = P_sig_watt + (abs(Dinf(i)))^2;
    i=i+1;
end

if Fsig2 > 0 % consider second tone only if present
 i= Nsig2 - floor(Nbin_sig/2); % signal power is calculated over Nbin +/- Nbin_sign/2
 while i <= (Nsig2 + floor(Nbin_sig/2))  
    P_sig_watt = P_sig_watt + (abs(Dinf(i)))^2;
    i=i+1;
 end
end

%% calculate power at DC
P_dc_watt = 0; % Power at DC
i=1;
while i <= 6  
    P_dc_watt = P_dc_watt + (abs(Dinf(i)))^2;
    i=i+1;
end

%% Total power in the BW
P_total_watt = 0; % total power in BW
i=Nstart;
while i <= Nend  
    P_total_watt = P_total_watt + (abs(Dinf(i)))^2;
    i=i+1;
end

P_noise_watt = P_total_watt - P_sig_watt;

% Add noise in signal bins
if Fsig2 > 0 % check if second signal tone is present
  avg_noise_per_bin = 2*P_noise_watt/(Nend-Nstart-2*Nbin_sig);%
  Noise_in_sigbin = 2*Nbin_sig*avg_noise_per_bin;
else
  avg_noise_per_bin = P_noise_watt/(Nend-Nstart-Nbin_sig);
  Noise_in_sigbin = Nbin_sig*avg_noise_per_bin;
end

P_noise_watt = P_noise_watt + Noise_in_sigbin;

SNDR = P_sig_watt/P_noise_watt;
SNDR_dB = 10*log10(SNDR);

%% calculate SFDR
% Remove signal power and calculate spectrum

Pfft_nosig = Dinf;
i= Nsig1 - floor(Nbin_sig/2); % signal power is calculated over Nbin +/- Nbin_sign/2
while i <= (Nsig1 + floor(Nbin_sig/2))
    Pfft_nosig(i) = sqrt(avg_noise_per_bin);
    i=i+1;
end

if Fsig2 > 0 % consider second tone only if present
 i= Nsig2 - floor(Nbin_sig/2); % signal power is calculated over Nbin +/- Nbin_sign/2
 while i <= (Nsig2 + floor(Nbin_sig/2))
    Pfft_nosig(i) = sqrt(avg_noise_per_bin);
    i=i+1;
 end
 
end

Pspur_max_inband = abs(max(Pfft_nosig(Nstart:Nend))); % in band spur magnitude
Pspur_max_OOBlo = abs(max(Pfft_nosig(6:Nstart-1))); % out of band spur. Lower side of BW
Pspur_max_OOBhi = abs(max(Pfft_nosig(Nend+1:N/2))); % out of band spur. Higher side of  BW
Pspur_max_OOB = max(Pspur_max_OOBlo,Pspur_max_OOBhi); % max out of band spur
P_sig_peak = max(abs(Dinf(Nsig1-Nbin_sig/2:Nsig1+Nbin_sig/2)));
SFDR_dB_inband = 20*log10(P_sig_peak/Pspur_max_inband);
SFDR_dB_OOB = 20*log10(P_sig_peak/Pspur_max_OOB);

%% calculate HD3 or IM3

if Fsig2 > 0
 Nhd3_lo = ceil(N*(2*Fsig2-Fsig1)/Fs); % bin number for HD3
 Nhd3_hi = ceil(N*(2*Fsig1-Fsig2)/Fs);
 HD3 = max(Pfft_nosig(Nhd3_lo),Pfft_nosig(Nhd3_hi));
 HD3_dBc = 20*log10(HD3/P_sig_peak);
else
 Nhd3_hi = ceil(N*(3*Fsig1)/Fs);
 HD3 = Pfft_nosig(Nhd3_hi);
 HD3_dBc = 20*log10(HD3/P_sig_peak);
end

end
