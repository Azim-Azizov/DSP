%Task3
clc; clear; close all;

fs = 2000;      
fc = 460;      
M = 20;          
N = M + 1;      
n = 0:M;        
ft = fc / fs;    

h_sinc = sinc(2 * ft * (n - M/2));


w_hamming = hamming(N)';


h_windowed = h_sinc .* w_hamming;

% Truncated and shifted sinc function
figure;
stem(n - M/2, h_sinc, 'filled', 'LineWidth', 1.5);
xlabel('n');
ylabel('Amplitude');
title('Truncated and Shifted Sinc Function');
grid on;

% Hamming Window
figure;
subplot(221);
stem(n - M/2, w_hamming, 'filled', 'LineWidth', 1.5);
xlabel('n');
ylabel('Amplitude');
title('Hamming Window');
grid on;

% Impulse Response of Hamming-Windowed Low-Pass Filter
subplot(222)
stem(n - M/2, h_windowed, 'filled', 'LineWidth', 1.5);
xlabel('n');
ylabel('Amplitude');
title('Impulse Response of Hamming-Windowed Low-Pass Filter');
grid on;


H_windowed = fft(h_windowed, 1024);
H_rect = fft(h_sinc, 1024);  % Compare with rectangular window
H_mag_windowed = abs(H_windowed);
H_mag_rect = abs(H_rect);
H_dB_windowed = 20 * log10(H_mag_windowed);
H_dB_rect = 20 * log10(H_mag_rect);
f = linspace(0, fs/2, length(H_windowed)/2);

% Frequency Response (Linear Scale)
subplot(223);
plot(f, H_mag_rect(1:end/2), '--', 'LineWidth', 1.5, 'DisplayName', 'Rectangular Window');
hold on;
plot(f, H_mag_windowed(1:end/2), 'LineWidth', 1.5, 'DisplayName', 'Hamming Window');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Response (Linear Scale)');
legend;
grid on;

% Frequency Response (dB Scale)
subplot(224)
plot(f, H_dB_rect(1:end/2), '--', 'LineWidth', 1.5, 'DisplayName', 'Rectangular Window');
hold on;
plot(f, H_dB_windowed(1:end/2), 'LineWidth', 1.5, 'DisplayName', 'Hamming Window');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Frequency Response (dB Scale)');
legend;
grid on;