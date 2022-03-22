function [DoutFFT] = FFT_compute(DinT)
% Returns FFT of Din after Blackmanharris window
N=length(DinT);
A_window = transpose(blackmanharris(N));
rms_w = sqrt(((norm(A_window,2))^2)/N);
mean_w = mean(A_window);
ECF = 1/rms_w;
ACF = 1/mean_w;
DinT_window = DinT.*A_window;
DoutFFT = abs(fft(DinT_window,N)*ACF*2/N);
end

