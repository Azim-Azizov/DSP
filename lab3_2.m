clc;
fs = 2000;      
fc = 460;        
M = 20;          
n = 0:M;      
ft = fc / fs;    
w_lpf= sin(2 * pi * ft * (n - M/2)) ./ (pi * (n - M/2));
w_lpf(M/2+1)=2*ft;

figure;
stem(n, w_lpf, 'b', 'LineWidth', 1.5); grid on;
title('Truncated shifted sinc function');
[H, w] = freqz(w_lpf, 1, 1024, fs);
figure;
subplot(211); plot(w, abs(H), 'r', 'LineWidth', 1.5); grid on;
title('Frequency response');
subplot(212); plot(w, 20*log10(abs(H)), 'r', 'LineWidth', 1.5); grid on;
title('Frequency response (dB)'); grid on;
