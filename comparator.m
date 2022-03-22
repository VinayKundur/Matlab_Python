function [ comp_out,Td ] = comparator( Vin_diff,Vos )
% SAR ADC Comparator Model
% Vin_diff = differential input voltage
% Vos = offset voltage

Vn = 100e-6 ; % input referred thermal Noise. 1-sigma
Td=0;

Vnoise = normrnd(0,Vn); % generate random number with mean and sigma  

if (Vin_diff+Vos+Vnoise)>0
    comp_out=1;
else
    comp_out=-1;
end


end

