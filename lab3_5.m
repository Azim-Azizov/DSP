%Task 5
clc; clear; close all;

fs = 2000;    
fc = 460;     
M = 20;       
N = M + 1;    
n = 0:M;      
ft = fc / fs; 

h_sinc = sinc(2 * ft * (n - M/2));
w_hamming = hamming(N)';
h_lp = h_sinc .* w_hamming;
h_hp = h_lp .* (-1).^n;

H_lp = fft(h_lp, 1024);
H_hp = fft(h_hp, 1024);

H_dB_lp = 20 * log10(abs(H_lp));
H_dB_hp = 20 * log10(abs(H_hp));

f = linspace(0, fs/2, length(H_lp)/2);

figure;
subplot(211);
stem(n - M/2, h_lp, 'filled', 'LineWidth', 1.5, 'Color', [0 0.5 1]);
title('Impulse Response of Low-Pass Filter');
xlabel('n'); ylabel('Amplitude'); grid on;

subplot(212);
stem(n - M/2, h_hp, 'filled', 'LineWidth', 1.5, 'Color', [1 0.2 0.2]);
title('Impulse Response of High-Pass Filter');
xlabel('n'); ylabel('Amplitude'); grid on;

figure;
subplot(211);
plot(f, abs(H_lp(1:end/2)), 'LineWidth', 1.5, 'Color', [0 0.5 1], 'DisplayName', 'Low-Pass Filter');
hold on;
plot(f, abs(H_hp(1:end/2)), '--', 'LineWidth', 1.5, 'Color', [1 0.2 0.2], 'DisplayName', 'High-Pass Filter');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Response (Linear Scale)');
legend; grid on; hold off;

subplot(212);
plot(f, H_dB_lp(1:end/2), 'LineWidth', 1.5, 'Color', [0 0.5 1], 'DisplayName', 'Low-Pass Filter');
hold on;
plot(f, H_dB_hp(1:end/2), '--', 'LineWidth', 1.5, 'Color', [1 0.2 0.2], 'DisplayName', 'High-Pass Filter');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Frequency Response (dB Scale)');
legend; grid on; hold off;