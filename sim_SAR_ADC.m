% simulate SAR ADC model of Drone Receiver baseband
% Author : Vinay Kundur
% Last Modified: 22nd Oct 2021
% verified in MATLAB 7.9.0 R2009b

tic;
clear;
clc;

%% Create results Directory
dir_path = 'C:\Users\kundurv\Documents\Projects-AIS\Drone\Drone_RX\MATLAB\SAR_ADC_model\Results\'; % results directory. Sub directory will be created for each run with time stamp. No need to change this for every run
dir_name = horzcat('droneRXADC_MATLAB_results_',datestr(clock,'yyyymmddHHMMSS'),'\'); % separate results folder for each run based on time stamp
results_dir = horzcat(dir_path,dir_name);
mkdir(results_dir);

logfile = horzcat(results_dir,'logfile.txt');
logfileID = fopen(logfile,'w');
fprintf(logfileID,'Simulation started on %s \n',datestr(clock));
fprintf(logfileID,'Simulation run by: %s \n', getenv('username'));

%% Set Primary variables for simulations
sim_description = 'SAR ADC simulations without calibration, with thermal noise'; % brief description sim condition to be added to log file
fprintf(logfileID,'Simulation description: %s \n', sim_description);

sim_iterations = 2; % number of sim iterations to cover cap and noise variations
Nfft = 2^12; % Number of FFT points
Fclk = 24.01e6; % sampling clock
Jph_rms = 50e-12; % sampling clock rms phase jitter
N = 10; % ADC resolution. SAR ADC model code needs to be changed if resolution is changed due to fixed DAC architecture.
fullscale = (2^N)-1; % Full scale of waveform. For 10 bit ADC, it is 1023
fullscale_V = 2; % full scale in volts. 2*VDD
commonmode = (2^N)/2; % common mode code of the signal that needs to be subtracted from output code
%%signal
Fin1 = 4.113e6 ; % 1st signal frequency. max = 4.4MHz
Vin1 = 2.5e-3 ; % 1st signal sine wave amplitude % -95dBm signal = 2.9mVrms
Fin2 = Fin1+251e3 ; % 2nd signal frequency. max = 4.4MHz
Vin2 = Vin1 ; % 2nd signal sine wave amplitude % -95dBm signal = 2.9mVrms
channel_BW = 1.5e6 ; % channel BW used for SNDR calculation. 1MHz
%%Blocker
Fblk1 = 6.01317e6 ; % blocker frequency
Vblk1 = 0e-3 ; % blocker power. -48dBm = 625mV
Fblk2 = 2*Fblk1-(Fin1+Fin2)/2 ; % 2nd blocker frequency chosen such that Fblk1,Fclk2 IMD fall inband
Vblk2 = 0e-3 ; % blocker power. -48dBm = 625mV
%%Non idealities
Vn_rms_prestg = 650e-6; % 641uV. RMS noise from LNA+MXR ( 580uV) and VGA+LPF(275uV)
Vn_comparator = 400e-6; % 400uV. SAR ADC comparator RMS thermal noise
Vos_comparator = 3.5e-3 ; % 3.5mV. comparator offset voltage 1-sigma
Vn_ref = 50e-6; % 50uV. reference voltage rms noise
cap_mismatch_1sigma = 0.5 ; % 1-sigma capacitor mismatch for 32fF cap in percent. 0.25% for 32fF from simulations
cal_en = 1; % enable 4 MSB cap calibration

SNDR_spec = 13; % SNDR specification at ADC output in dB

%% update derived variables
Fstart= (Fin1 + Fin2)/2 - channel_BW/2; % Bandwidth start frequency
Fend = (Fin1 + Fin2)/2 + channel_BW/2; % Bandwidth end frequency

pi = 3.1415926535;
t=0:1/Fclk:(Nfft-1)/Fclk;
t_err = normrnd(0,Jph_rms,[1,length(t)]);
t=t+t_err;
sig1 = Vin1*sin(2*pi*Fin1*t);
sig2 = Vin2*sin(2*pi*Fin2*t);
blk1 = Vblk1*sin(2*pi*Fblk1*t);
blk2 = Vblk2*sin(2*pi*Fblk2*t);

Vn = normrnd(0,Vn_rms_prestg,[1,length(t)]); % noise from LNA+MIXER+VGA+LPF at ADC input

sig = sig1 + sig2 + blk1 + blk2 + Vn; % signal at ADC input

%% plot input time domain waveform

figure(1);
subplot(2,1,1);
plot(t*1e9,sig,'LineWidth',3);
set(gca,'FontSize',12);
grid on;
ylabel('input signal(V)');
xlabel('Time(ns)');
set(gcf,'color','w');

subplot(2,1,2); % zoomed in version of time domain plot
Nplot = 40; % number of samples for plotting
plot(t(10:10+Nplot)*1e9,sig(10:10+Nplot),'LineWidth',3);
set(gca,'FontSize',12);
grid on;
legend('Zoomed in waveform');
ylabel('input signal(V)');
xlabel('Time(ns)');
set(gcf,'color','w');
saveas(gcf,horzcat(results_dir,'input_waveform.png'));

%% Signal goes through SAR ADC model
%caps_all = zeros(sim_iterations,N);
DNL_all = zeros(sim_iterations,2^N);
INL_all = zeros(sim_iterations,2^N);
SNDR_inband_all = zeros(1,sim_iterations);
SFDR_inband_all = zeros(1,sim_iterations);
SFDR_OOB_all = zeros(1,sim_iterations);
Dout_all = zeros(sim_iterations,Nfft);

A_window = transpose(blackmanharris(Nfft));
W1 = norm(A_window,1);
W2 = norm(A_window,2);
NBW = (W2/W1)^2;
NRBW = NBW*Fclk;
rms_w = sqrt(((norm(A_window,2))^2)/Nfft);
mean_w = mean(A_window);
ECF = 1/rms_w;
ACF = 1/mean_w; % DC value of the window
CF = ECF;

fprintf(logfileID,'\nTotal Iterations:%d \nIteration progress:\n',sim_iterations);
fprintf('Total Iterations:%d \nIteration progress:\n',sim_iterations);

for i=1:1:sim_iterations
    [Dout_all(i,:),caps_all(i,:),DNL_all(i,:),INL_all(i,:),code_ramp(i,:),Vin_RAMP,caps_caldac(i,:),cal_code(i,:),VFS_sim(i),Vos_comp(i)] = SAR_ADC_redundancy_2BridgeCAP_digCAL(sig,N,Fclk,cap_mismatch_1sigma,cal_en,Vos_comparator,Vn_comparator,Vn_ref,fullscale_V/2,results_dir,i); % signal going through SAR ADC
    DinT_window = (Dout_all(i,:)-commonmode).*A_window;
    DoutFFT = abs(fft(DinT_window,Nfft)*CF*4/Nfft); % FFT/N returns A/2 ( two sided spectrum). For single sided,*2 and for a fullscale sine FS, amplitude is FS/2, so *2 for dBFS
    [SNDR_inband_all(i),SFDR_inband_all(i),SFDR_OOB_all(i),~,~,~ ,~]= SNDR(DoutFFT,Fin1,Fin2, Fclk, Fstart,Fend);
    fprintf(logfileID,'%d out of %d\n',i,sim_iterations);
    fprintf('%d out of %d\n',i,sim_iterations);
    
end
fprintf('\n');

num_fails = sum(SNDR_inband_all < SNDR_spec); % number of samples with SNDR < 19.5;
pass_percent = (sim_iterations-num_fails)*100/sim_iterations;

Dout = Dout_all(find(SNDR_inband_all == min(SNDR_inband_all)),:); % copy data of the iterations where SNDR is minimum

figure(2); % plotting DNL/INL
subplot(2,1,1);
plot(0:1:(2^N)-1,DNL_all,'-','LineWidth',3);
set(gca,'FontSize',12);
set(gcf,'color','w');
grid on;
xlim([0 (2^N)-1]);
% ylim([-2 2]);
ylabel('DNL(LSB)');
xlabel('CODE');

subplot(2,1,2);
plot(0:1:(2^N)-1,INL_all,'-','LineWidth',3);
set(gca,'FontSize',12);
set(gcf,'color','w');
grid on;
xlim([0 (2^N)-1]);
% ylim([-2 2]);
ylabel('INL(LSB)');
xlabel('CODE');
saveas(gcf,horzcat(results_dir,'DNL_INL.png'));

figure(3); % Plot input voltage Vs output code
plot(Vin_RAMP,code_ramp,'-','LineWidth',3);
set(gca,'FontSize',12);
set(gcf,'color','w');
grid on;
% xlim([-1 1]);
% ylim([0 (2^N)-1]);
ylabel('ADC output code');
xlabel('Input Voltage(V)');
saveas(gcf,horzcat(results_dir,'input_vs_output.png'));

SNDR_mean = mean(SNDR_inband_all);
%SNDR_min = min(SNDR_inband_all);
%SNDR_max = max(SNDR_inband_all);
SNDR_sigma = std(SNDR_inband_all);

SFDR_OOB_mean = mean(SFDR_OOB_all);
%SFDR_OOB__min = min(SFDR_OOB_all);
%SFDR_OOB__max = max(SFDR_OOB_all);
SFDR_OOB_sigma = std(SFDR_OOB_all);

%% plot output FFT
DinT_window = (Dout-commonmode).*A_window;
DoutFFT = abs(fft(DinT_window,Nfft)*CF*4/Nfft); % FFT/N returns A/2 ( two sided spectrum). For single sided,*2 and for a fullscale sine FS, amplitude is FS/2, so *2 for dBFS
fullscale_dB = 20*log10(fullscale);
DoutFFT_dB = 20*log10(DoutFFT)-fullscale_dB; % For normalizing w.r.t Fullscale
frequency = Fclk/Nfft:Fclk/Nfft:Fclk;

[SNDR_dB,SFDR_inband_dB,SFDR_OOB_dB,HD3_dBc,Psig_watt,Pnoise_watt ,DoutFFTnosig]= SNDR(DoutFFT,Fin1,Fin2, Fclk, Fstart,Fend);
Vi_offset_mV = (mean(Dout)-commonmode)*VFS_sim(end)*1000/fullscale;

if Fin2 > 0
    sig_rms_LSB = sqrt(Psig_watt)*0.5;
    noise_LSB = sqrt(Pnoise_watt)*0.5;
    Vsig_rms_mV = sqrt(Psig_watt)*1e3*0.5*fullscale_V/(sqrt(2)*fullscale);
    Vn_rms_uV = sqrt(Pnoise_watt)*1e6*0.5*fullscale_V/(sqrt(2)*fullscale);
else
    sig_rms_LSB = sqrt(Psig_watt)*0.707;
    noise_LSB = sqrt(Pnoise_watt)*0.707;
    Vsig_rms_mV = sqrt(Psig_watt)*1e3*0.707*fullscale_V/(sqrt(2)*fullscale);
    Vn_rms_uV = sqrt(Pnoise_watt)*1e6*0.707*fullscale_V/(sqrt(2)*fullscale);
end

DoutFFTnosig_dB = 20*log10(DoutFFTnosig)-fullscale_dB; % For normalizing w.r.t Fullscale

figure(4); % DFT
plot(frequency(1:Nfft/2)*1e-6,DoutFFT_dB(1:Nfft/2),'r','LineWidth',3);
hold on;
plot(frequency(1:Nfft/2)*1e-6,DoutFFTnosig_dB(1:Nfft/2),'b','LineWidth',3);
hold off;
title('Worst Case SNDR corner');
legend('with signal','signal removed');
set(gca,'FontSize',12);
grid on;
ylabel('PSD(dBFS/NBW)');
ylim([-100 0]);
textdisplay = horzcat('(RBW=', num2str(NRBW*1e-3), 'KHz)');
xlabel(horzcat('Frequency(MHz) ',textdisplay));
textdisplay = horzcat('SNDR=', num2str(min(SNDR_inband_all)), 'dB');
text(2,-20,textdisplay,'FontSize',14);
xlim([0 Fclk*(1e-6)/2]);
set(gcf,'color','w');
saveas(gcf,horzcat(results_dir,'FFT_worstcase.png'));

figure(5); % plot SNDR Vs number of samples
spec = SNDR_spec*ones(1,sim_iterations); % spec is 20dB
plot(1:1:sim_iterations, SNDR_inband_all,1:1:sim_iterations, spec,'-','LineWidth',3);
set(gca,'FontSize',12);
set(gcf,'color','w');
grid on;
% xlim([-1 1]);
% ylim([0 (2^N)-1]);
ylabel('SNDR Inband(dB)');
xlabel('Sample Number');
legend('SNDR across samples','SNDR spec');
textdisplay = horzcat('Yield=', num2str(pass_percent), '%');
text(1,SNDR_spec+2,textdisplay,'FontSize',14);
textdisplay = horzcat('Mean=', num2str(mean(SNDR_inband_all)), 'dB');
text(1,SNDR_spec+1.5,textdisplay,'FontSize',14);
textdisplay = horzcat('std dev=', num2str(std(SNDR_inband_all)), 'dB');
text(1,SNDR_spec+1,textdisplay,'FontSize',14);
saveas(gcf,horzcat(results_dir,'SNDR_vs_samples.png'));

figure(6); % plot histogram of SNDR
spec = SNDR_spec*ones(1,sim_iterations); 
hist(SNDR_inband_all);
set(gca,'FontSize',12);
set(gcf,'color','w');
grid on;
xlim([10 25]);
% ylim([0 (2^N)-1]);
ylabel('# of samples');
xlabel('SNDR inband(dB)');
textdisplay = horzcat('Yield=', num2str(pass_percent), '%');
text(SNDR_spec,sim_iterations/5,textdisplay,'FontSize',14);
textdisplay = horzcat('Mean=', num2str(mean(SNDR_inband_all)), 'dB');
text(SNDR_spec,sim_iterations/6,textdisplay,'FontSize',14);
textdisplay = horzcat('std dev=', num2str(std(SNDR_inband_all)), 'dB');
text(SNDR_spec,sim_iterations/7,textdisplay,'FontSize',14);
line([SNDR_spec, SNDR_spec],ylim,'LineWidth',2,'Color','r');
saveas(gcf,horzcat(results_dir,'SNDR_histogram.png'));



figure(7); % plot calibration codes
plot(1:1:sim_iterations, cal_code,'-','LineWidth',3);
set(gca,'FontSize',12);
set(gcf,'color','w');
grid on;
% xlim([-1 1]);
% ylim([0 255]); % CALCODE is 8 bit
ylabel('Calibration Code');
xlabel('Sample Number');
legend('MSB1 CALCODE','MSB2 CALCODE','MSB3 CALCODE','MSB4 CALCODE','MSB5 CALCODE');
saveas(gcf,horzcat(results_dir,'CALCODES.png'));


fprintf(logfileID,'\n\nSimulations ended on %s\n',datestr(clock));
save(horzcat(results_dir,'variables.mat'));
fprintf(logfileID,'\nTotal simulation time = %f seconds\n',toc);
fclose('all');
toc;
