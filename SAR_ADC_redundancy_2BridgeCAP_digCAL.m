function [Dout,caps,DNL,INL,code_ramp,Vin_RAMP,caps_caldac,rel_cal_code_corrected,VFS_sim,Vos_comp] = SAR_ADC_redundancy_2BridgeCAP_digCAL(Vin,N,Fclk,Cmis_sigma_percent,EN_CAL,Vos_comp_sigma,Vn_comparator, Vref_noise_rms, VR, results_dir,sim_num)
%% DAC with redundancy and digital calibration.
% Based on ISSCC 2010 21.5 C.C Liu "A 10b 100Msps 1.13mW SAR ADC with
% binary scaled Error compensation"
% But location of redundant weights are changed
% 512 256 128
%         128 64 32
%                32 16
%                   16 8 4 2 1 0.5 0.25 0.25
% Bridge1 cap added between 128 and 64 for cap scaling
% Bridge2 cap added between   8 and  4 for cap scaling
% 1 is the LSB
% 0.5 and 0.25 are present for calibration


%% Model Inputs
% Vin = input signal
% N = number of ADC bits (10)
% Fclk = sampling frequency (24MHz)
% Cmis_sigma_percent = 32fF cap mismatch 1-sigma in %
% EN_CAL = Enable calibration. 0 or 1
% Vos_comp_sigma = comaprator offset 1-sigma ( 3.5mV)
% Vref_noise_rms = Reference voltage noise rms
% VR = reference voltage value 1V
% results_dir = results directory
% sim_num = simulation iteration number

Rsw = 200 ; % sampling cap switch resistance
K = 1.38064852e-23 ; % Bolzmann constant
T = 300 ; % Temperature in kelvin

R = 3 ;% number of reduntant states
Ncc1 = 4; % location of compensation caps
Ncc2 = 7;
Ncc3 = 9;
Nr = N + R; % effective raw bits from ADC before digital correction
wt_red = 178; % total redundant weights

filename = horzcat(results_dir,'data_sample',num2str(sim_num),'.csv');
fileID = fopen(filename,'w');

%% Assign Cap Values and add mismatch

% binary equivalent array = [512 256 128 128 64 32 32 16 16 8 4 2 1 0.5 0.25 0.25]

% configuration 1
scale_cap = 160/512; % scale all the capacitors according to the need. 160/512 means MSB cap is 160fF

caps_upper = scale_cap*[512 256 128 128].*1e-15; % caps on MSB side of bridge cap

cbridge1 = 1*scale_cap*202.2857143e-15;
cpex1 = 0e-15; % pex cap across cbridge1 cap. This will create imbalance in cap wieghts
caps_mid = scale_cap*[512 256 256 128 128 64].*1e-15; % caps on LSB side of bridge cap

cbridge2 = scale_cap*76.8e-15; % bridge cap between mid and lower cap array.Cbridge = Ctotal(lower)*Ctotal(binary equivalent)/(Ctotal(lower)-Ctotal(binary equivalent))
cpex2 = 0e-15; % pex cap across cbridge1 cap. This will create imbalance in cap wieghts
ctrim = scale_cap*128; % trim cap to adjust error in cbridge2
caps_lower = scale_cap*[512 256 128 64 32 32 ctrim/scale_cap].*1e-15; % caps on LSB side of bridge cap


% apply cap mismatch
for i=1:1:length(caps_upper) % apply mismatch to caps on MSB side of bridge cap
    sigma_percent = Cmis_sigma_percent/sqrt(caps_upper(i)/32e-15); % scale cap mismatch sigma according to its size w.r.t 32fF cap
    caps_upper(i) = caps_upper(i)*(1 + normrnd(0,sigma_percent/100));
end

sigma_percent = Cmis_sigma_percent/sqrt(cbridge1/32e-15); % scale cap mismatch sigma according to its size w.r.t 32fF cap
cbridge1 = cbridge1*(1 + normrnd(0,sigma_percent/100));

for i=1:1:length(caps_mid)
    sigma_percent = Cmis_sigma_percent/sqrt(caps_mid(i)/32e-15); % scale cap mismatch sigma according to its size w.r.t 32fF cap
    caps_mid(i) = caps_mid(i)*(1 + normrnd(0,sigma_percent/100));
end

sigma_percent = Cmis_sigma_percent/sqrt(cbridge2/32e-15); % scale cap mismatch sigma according to its size w.r.t 32fF cap
cbridge2 = cbridge2*(1 + normrnd(0,sigma_percent/100));

for i=1:1:length(caps_lower)
    sigma_percent = Cmis_sigma_percent/sqrt(caps_lower(i)/32e-15); % scale cap mismatch sigma according to its size w.r.t 32fF cap
    caps_lower(i) = caps_lower(i)*(1 + normrnd(0,sigma_percent/100));
end

cbridge1 = cbridge1 + cpex1;
cbridge2 = cbridge2 + cpex2;


cap_scale_factor2 = cbridge2/(cbridge2+sum(caps_lower)); % cap scaling due to Cbridge2
caps_mid_low = [caps_mid caps_lower*cap_scale_factor2]; % combine mid and lower caps using cbridge2 scaling factor.
cap_scale_factor1 = cbridge1/(cbridge1+sum(caps_mid_low)); % cap scaling due to Cbridge1
caps = [caps_upper caps_mid_low*cap_scale_factor1]; % combine all caps into equivalent binary weights withtout bridge caps

% %from schematic
% caps = [174.4 87.2  43.6  43.6  21.9  10.95 10.95 5.474 5.474 2.737 1.393 0.696 0.348 0.174 0.086 0.086].*1e-15;
% %from PEX
% caps = [179.2 91.68 47.25 46.94 24.58 13.63 14.11 8.43  8.4   5.456 4.009 3.3   2.97  2.803 2.719 2.719].*1e-15;

fprintf(fileID, 'DAC capacitors (equivalent caps without bridge caps)\n');
fprintf(fileID,'%.6e,', caps);
%% Comparator offset
Vos_comp = normrnd(0,Vos_comp_sigma); % comparator offset
% Vos_comp = 18e-3;

ctotal = sum(caps);
Csample = ctotal; % Sampling cap. Total cap in DAC.

%% CAP mismatch calibration routine
% binary equivalent array = [512 256 128 128 64 32 32 16 16 8 4 2 1 0.5 0.25 0.25 ctrim]
N_msb_cal = 4; % number of MSB to calibrate. Max = 6.
N_cal_avg = 8; % number of times calcode is averaged to get final value. Min=2 and then powers of 2 for higher value.
caps_caldac = caps(end-8:end-1); % 5 LSB and 3 sub LSB of main DAC are used as CALDAC
N_caldac = length(caps_caldac);


Cref(1) = caps(2)+caps(3)+caps(4); % reference cap for calibrating MSB1
Cref(2) = caps(3)+caps(4); % reference cap for calibrating MSB2
Cref(3) = caps(4); % reference cap for calibrating MSB3
Cref(4) = caps(5)+caps(6)+ caps(7); % reference cap for calibrating MSB4
Cref(5) = caps(6)+ caps(7); % reference cap for calibrating MSB5
Cref(6) = caps(7); % reference cap for calibrating MSB6

rel_cal_code_corrected = zeros(1,6);
rel_cal_code = zeros(1,6);
LUT = zeros(1,2^N_msb_cal);

if EN_CAL == 1
    
    fprintf(fileID, '\n\n8 bit data from each iteration ( bits inverted in step 0 2 4 6 to account for VREF inversion)');
    for i=1:1:N_msb_cal % calibrate each MSB one by one
        C_cal = caps(i); % cap to be calibrated
        C_ref = Cref(i); % reference cap for calibration
        bits = zeros(1,N_caldac); % 1-sign bit.5 magnitude bits
        for x=1:1:N_cal_avg/2
            for step=1:1:2
                % step1:connect C_cal to VREFP and C_ref to VREFM
                % step2:connect C_cal to VREFM and C_ref to VREFP
                if step==1
                    vin_comp = ( C_cal*VR - C_ref*VR)/ctotal;
                else
                    vin_comp = (-C_cal*VR + C_ref*VR)/ctotal;
                end
                
                for j=1:1:N_caldac
                    [decision,Td]=comparator(vin_comp,Vos_comp); % comparator decision is +/-1
                    
                    if(decision>0)
                        bits(j)=1;
                    else
                        bits(j)=0;
                    end
                    
                    
                    if j <= N_caldac
                        vin_comp = vin_comp - (decision*VR*caps_caldac(j)/ctotal)+ normrnd(0,Vref_noise_rms) + normrnd(0,sqrt(K*T/Csample)); % switch caldac caps one by one
                    end
                end
                
                if(step==1)  % invert the bits of step1
                    for j=1:1:length(bits)   % bit =1 connect to +VREF and bit=0 means connect to -VREF
                        if bits(j)== 0
                            bits(j)= 1;
                        else
                            bits(j)= 0;
                        end
                    end
                end
                
                
                code(step,x) = bin2dec(strcat(num2str(bits(1:end)))); % weighted error code
            end
            
        end
        fprintf(fileID, '\n\nMSB Number %d \n',i);
        fprintf(fileID,'%d,', code);
        
        cal_code_raw(i) = nearest(sum(sum(code))/N_cal_avg); % calibration code is the average of step1 and step2 and number of iterations
        fprintf(fileID, '\nMSB Number %d code average\n',i);
        fprintf(fileID,'%d,', cal_code_raw(i));
    end
    
    
    % cal code =31 and 32 both are not adding any correction cap. + side and
    % -side connected caps cancel out
    %  effective magnitude is cal_code-132.
    

    
    for i =1:1:N_msb_cal
        if cal_code_raw(i) < 128
            rel_cal_code(i) = cal_code_raw(i) - 127; % relative cal code w.r.t Mid code
        else
            rel_cal_code(i) = cal_code_raw(i) - 128; % relative cal code w.r.t Mid code
        end
    end
    
    fprintf(fileID, '\n\nRelative calcode for MSBs:\n');
    fprintf(fileID,'%d,', rel_cal_code);
    
    % binary equivalent array = [512 256 128 128 64 32 32 16 16 8 4 2 1 0.5 0.25 0.25]
    % applying the error code of the reference cap value
    rel_cal_code_corrected(6) = sum(rel_cal_code(6));
    rel_cal_code_corrected(5) = sum(rel_cal_code(5)+rel_cal_code_corrected(6));
    rel_cal_code_corrected(4) = sum(rel_cal_code(4)+rel_cal_code_corrected(5)+rel_cal_code_corrected(6));
    rel_cal_code_corrected(3) = sum(rel_cal_code(3)+rel_cal_code_corrected(4));
    rel_cal_code_corrected(2) = sum(rel_cal_code(2)+rel_cal_code_corrected(3)+ rel_cal_code_corrected(4));
    rel_cal_code_corrected(1) = sum(rel_cal_code(1)+rel_cal_code_corrected(2)+ rel_cal_code_corrected(3)+ rel_cal_code_corrected(4));
    
    fprintf(fileID, '\n\nRelative calcode for 4 MSB after correcting for errors in reference caps:\n');
    fprintf(fileID,'%d,', rel_cal_code_corrected);
    
    fprintf(fileID, '\n\nLook Up Table(LUT) Entries:\n');
    for i=1:1:length(LUT)
     bin_i = dec2bin(i-1,N_msb_cal);
     LUT_entry = 0;
        for j=1:1:length(bin_i)
         LUT_entry = LUT_entry + rel_cal_code_corrected(j)*str2num(bin_i(j));
        end
     LUT_entry = floor(LUT_entry/4 + 0.5); % CALDAC is 1/4 of main LSB. So, LUT entry corrected to match main DAC LSB
     LUT(i) = LUT_entry;
     fprintf(fileID,'\nLUT Addr %d,%d',i-1, LUT_entry);
    end
    
       
end

%% Emulate track and hold. Add KT/C noise
T_TH = 0.5/Fclk ; % input signal tracking time 0.5/Fs.
Vin_SH = zeros(1,length(Vin));
Vin_SH(1)=Vin(1);
Vn_sigma = sqrt(K*T/Csample);% KT/C noise Vrms

for i=2:1:length(Vin)
    Vin_SH(i) = Vin(i-1) + (Vin(i)-Vin(i-1))*(1-2.718281^(-T_TH/(Rsw*Csample))); % Model RC settling for track and hold
    Vin_SH(i) = Vin_SH(i) + normrnd(0,Vn_sigma); % Add KT/C noise
end


%% Compute DNL and INL through ramp test
Vstep = 2*VR/((2^N-1)+wt_red);
Vin_RAMP = -VR:Vstep:VR;
%Vin_RAMP = Vin_RAMP - Vstep/4 ; % offset input by LSB/2 so that decision threshold and LSB are off by LSB/2
code_ramp = zeros(1,length(Vin_RAMP));


vin_comp = 0;
fprintf(fileID,'\n\nInput Voltage(V), 13bit raw code(dec), 10bit final code(dec)');
for i=1:1:length(Vin_RAMP)  % for each voltage sample
    vin_comp = Vin_RAMP(i);
    bits = zeros(1,Nr);
    correction_code = 0;
    for j=1:1:Nr
        [decision,Td]=comparator(vin_comp,Vos_comp); % comparator decision is +/1
        if Td > (0.5/(Fclk*Nr))  % if comparator delay is high due to metastability.Just keep output to zero.
            decision=-1;
        end
        
        Ceff = caps(j)*(ctotal-caps(j))/(caps(j)*(ctotal-caps(j))); % effective cap for noise calculation. New switched cap is in series with rest of the caps
        noise_gain =  caps(j)/ctotal; % noise voltage division due to capacitors
        Vn_KTC_SARstep = noise_gain*sqrt(K*T/Ceff);
        
        vin_comp = vin_comp - (decision*VR*(caps(j)/ctotal))+ normrnd(0,Vref_noise_rms) + normrnd(0,Vn_KTC_SARstep); % scale comparator input based on decision and add VREF noise and KT/C noise
        
        
        if(decision>0)
            bits(j)=1;
        else
            bits(j)=0;
        end        
    end
    
    %convert raw redundant bits to binary outputs
    part_lsb = strcat(num2str(bits(Ncc3:end)));
    part_lo = strcat(num2str(bits(Ncc2:Ncc3-1)));
    part_mid = strcat(num2str(bits(Ncc1:Ncc2-1)));
    part_hi = strcat(num2str(bits(1:Ncc1-1)));
    correction_code = LUT(bin2dec(strcat(num2str(bits(1:4)))) + 1);
    code_ramp(i) = bin2dec(part_hi)*(2^(N-Ncc1+1))+ bin2dec(part_mid)*(2^(N-Ncc2+2)) + bin2dec(part_lo)*2^(N-Ncc3+3)+ bin2dec(part_lsb)-(2^(N-Ncc1) + 2^(N-Ncc2+1)+ 2^(N-Ncc3+2)) - correction_code;
    
    
        if code_ramp(i) > 1023
            code_ramp(i) = 1023;
        elseif code_ramp(i) < 0
            code_ramp(i) = 0;
        else
            code_ramp(i) = code_ramp(i);
        end
    fprintf(fileID,'\n%f, %d, %d',Vin_RAMP(i), bin2dec(strcat(num2str(bits(1:end)))),code_ramp(i));
    
end
LSB_step = zeros(1,2^N);
LSB_step(2:end) = code_ramp(0.5*wt_red+2:end-0.5*wt_red)-code_ramp(0.5*wt_red+1:end-0.5*wt_red-1);
LSB_step(1)=mean(LSB_step(2:end));
%LSB_step(1)= 1;
DNL = (LSB_step - LSB_step(1))/LSB_step(1);
INL = cumsum(DNL);
vin_code0 = Vin_RAMP(find(code_ramp>0, 1 ));
vin_code1023 = Vin_RAMP(find(code_ramp<1024, 1, 'last' ));
VFS_sim = vin_code1023 - vin_code0;


%% SAR Conversions of input signal
Dout = zeros(1,length(Vin));
vin_comp = 0;
for i=1:1:length(Vin_SH)  % for each voltage sample
    vin_comp = Vin_SH(i);
    bits = zeros(1,Nr);
    correction_code = 0;
    for j=1:1:Nr
        [decision,Td]=comparator(vin_comp,Vos_comp); % comparator decision is +/1
        if Td > (0.5/(Fclk*Nr))  % if comparator delay is high due to metastability.Just keep output to zero.
            decision=-1;
        end
        
        Ceff = caps(j)*(ctotal-caps(j))/(caps(j)*(ctotal-caps(j))); % effective cap for noise calculation. New switched cap is in series with rest of the caps
        noise_gain =  caps(j)/ctotal; % noise voltage division due to capacitors
        Vn_KTC_SARstep = noise_gain*sqrt(K*T/Ceff);
        
        
        vin_comp = vin_comp - (decision*VR*(caps(j)/ctotal))+ normrnd(0,Vref_noise_rms) + normrnd(0,Vn_KTC_SARstep); % scale comparator input based on decision and add VREF noise and KT/C noise
        
        
        if(decision>0)
            bits(j)=1;
        else
            bits(j)=0;
        end
        
        if j <= N_msb_cal
            if(decision > 0)
                correction_code = correction_code + rel_cal_code_corrected(j); % apply correction code only if bit=1
            end
        end
    end
    
    part_lsb = strcat(num2str(bits(Ncc3:end)));
    part_lo = strcat(num2str(bits(Ncc2:Ncc3-1)));
    part_mid = strcat(num2str(bits(Ncc1:Ncc2-1)));
    part_hi = strcat(num2str(bits(1:Ncc1-1)));
    correction_code = LUT(bin2dec(strcat(num2str(bits(1:4)))) + 1);
    Dout(i) = bin2dec(part_hi)*(2^(N-Ncc1+1))+ bin2dec(part_mid)*(2^(N-Ncc2+2)) + bin2dec(part_lo)*2^(N-Ncc3+3)+ bin2dec(part_lsb)-(2^(N-Ncc1) + 2^(N-Ncc2+1)+ 2^(N-Ncc3+2)) - correction_code;
    
        if Dout(i) > 1023
            Dout(i) = 1023;
        elseif Dout(i) < 0
            Dout(i) = 0;
        else
            Dout(i) = Dout(i);
        end
    
end

fclose(fileID);

end
