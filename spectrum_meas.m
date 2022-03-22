%Spectrum measurement
%Vinay Kundur
%Last updated: 28th March 2021

clear; clc;
fullscale = 1023; % Full scale of waveform. For 10 bit ADC, it is 1023
fullscale_V = 2; % full scale in volts
commonmode = 512; % common mode of the signal that needs to be subtracted from output
Fsig1=4.4123e6; % 1st tone signal frequency. 4.4123MHz
Fsig2=3.5103e6; % 2nd tone signal frequency. Set this to 0 for single tone. 3.5103MHz
Fstart=3.41e6; % Bandwidth start frequency
Fend=4.51e6; % Bandwidth end frequency

%% read CSV file
filepath = 'C:\Users\kundurv\Documents\Projects-AIS\Drone\Drone_RX\MATLAB\CSV_files\Async\';
filename = 'ADCLPF_m95dBm_Pblk55dBm6MHz_Therm5G_Jph25ps_TTSSFF_pexOTACOMP.csv';
data = csvread(horzcat(filepath,filename),1,0);
time = data(:,1);
voltage = data(:,4);
voltage=voltage-commonmode;
Tsample = time(2)-time(1);
Fsample = 1/Tsample;


%% plot time domain waveform
figure(1);
subplot(2,1,1);
plot(time*1e9,voltage,'LineWidth',3);
set(gca,'FontSize',12);
grid on;
ylabel('Voltage(V)');
xlabel('Time(ns)');
set(gcf,'color','w');


subplot(2,1,2); % zoomed in version of time domain plot
Nplot = 20; % number of samples for plotting
plot(time(10:10+Nplot,1)*1e9,voltage(10:10+Nplot,1),'LineWidth',3);
set(gca,'FontSize',12);
grid on;
legend('Zoomed in waveform');
ylabel('Voltage(V)');
xlabel('Time(ns)');
set(gcf,'color','w');


%% plot FFT
N=length(voltage);
A_window = blackmanharris(N);
W1 = norm(A_window,1);
W2 = norm(A_window,2);
NBW = (W2/W1)^2;
NRBW = NBW*Fsample;
rms_w = sqrt(((norm(A_window,2))^2)/N);
mean_w = mean(A_window);
ECF = 1/rms_w;
ACF = 1/mean_w; % DC value of the window
DinT_window = voltage.*A_window;
DoutFFT = abs(fft(DinT_window,N)*ECF*4/N); % FFT/N returns A/2 ( two sided spectrum). For single sided,*2 and for a fullscale sine FS, amplitude is FS/2, so *2 for dBFS
fullscale_dB = 20*log10(fullscale);
DoutFFT_dB = 20*log10(DoutFFT)-fullscale_dB; % For normalizing w.r.t Fullscale
frequency = Fsample/N:Fsample/N:Fsample;

[SNDR,SFDR_dB_inband,SFDR_dB_OOB,HD3_dBc,Psig_watt,Pnoise_watt ,DoutFFTnosig]= SNDR(DoutFFT,Fsig1,Fsig2, Fsample, Fstart,Fend);
Vi_offset_mV = mean(voltage)*fullscale_V*1000/fullscale;



if Fsig2 > 0
  sig_rms_LSB = sqrt(Psig_watt)*0.5;
  noise_LSB = sqrt(Pnoise_watt)*0.5;
  Vsig_rms_mV = sqrt(Psig_watt)*1e3*0.5*fullscale_V/fullscale;
  Vn_rms_uV = sqrt(Pnoise_watt)*1e6*0.5*fullscale_V/fullscale;
else
  sig_rms_LSB = sqrt(Psig_watt)*0.707;
  noise_LSB = sqrt(Pnoise_watt)*0.707;
  Vsig_rms_mV = sqrt(Psig_watt)*1e3*0.707*fullscale_V/fullscale;
  Vn_rms_uV = sqrt(Pnoise_watt)*1e6*0.707*fullscale_V/fullscale;
end

DoutFFTnosig_dB = 20*log10(DoutFFTnosig)-fullscale_dB; % For normalizing w.r.t Fullscale

figure(2);
plot(frequency(1:N/2)*1e-6,DoutFFT_dB(1:N/2),'r','LineWidth',3);
hold on;
plot(frequency(1:N/2)*1e-6,DoutFFTnosig_dB(1:N/2),'b','LineWidth',3);
hold off;
legend('with signal','signal removed');
set(gca,'FontSize',12);
grid on;
ylabel('PSD(dBFS/NBW)');
ylim([-100 0]);
textdisplay = horzcat('(RBW=', num2str(NRBW*1e-3), 'KHz)');
xlabel(horzcat('Frequency(MHz) ',textdisplay));
textdisplay = horzcat('SNDR(dB) inband=', num2str(SNDR), 'dB');
text(1,-20,textdisplay,'FontSize',14);
textdisplay = horzcat('SFDR(dB) inband=', num2str(SFDR_dB_inband), 'dB');
text(1,-30,textdisplay,'FontSize',14);
xlim([0 Fsample*(1e-6)/2]);
set(gcf,'color','w');


