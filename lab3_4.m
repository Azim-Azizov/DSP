%task4
clc; clear; close all;
fs = 44100;
fc = 5000;
M = 28; % Filter order
N = M + 1;
ft = fc / fs;
n = -M/2:M/2;
h = sinc(2 * ft * n); % Truncated sinc function
% Different Windows
w_rect = rectwin(N)';
w_ham = hamming(N)';
w_hann = hann(N)';
w_black = blackman(N)';
h_rect = h .* w_rect;
h_ham = h .* w_ham;
h_hann = h .* w_hann;
h_black = h .* w_black;
figure;
subplot(2,2,1); stem(n, w_rect, 'm'); title('Rectangular Window');
subplot(2,2,2); stem(n, w_ham, 'c'); title('Hamming Window');
subplot(2,2,3); stem(n, w_hann, 'y'); title('Hann Window');
subplot(2,2,4); stem(n, w_black, 'b'); title('Blackman Window');
% Frequency Responses
figure;
subplot(211);
hold on;
[H_rect, f] = freqz(h_rect, 1, 1024, fs);
[H_ham, ~] = freqz(h_ham, 1, 1024, fs);
[H_hann, ~] = freqz(h_hann, 1, 1024, fs);
[H_black, ~] = freqz(h_black, 1, 1024, fs);
plot(f, 20*log10(abs(H_rect)), 'm', 'DisplayName', 'Rectangular Window');
plot(f, 20*log10(abs(H_ham)), 'c', 'DisplayName', 'Hamming Window');
plot(f, 20*log10(abs(H_hann)), 'y', 'DisplayName', 'Hann Window');
plot(f, 20*log10(abs(H_black)), 'b', 'DisplayName', 'Blackman Window');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Comparison of Low-Pass Filters Using Different Windows');
legend;
grid on;
hold off;
% Plot Frequency Response (Non-dB)
subplot(212);
hold on;
plot(f, abs(H_rect), 'm', 'DisplayName', 'Rectangular Window');
plot(f, abs(H_ham), 'c', 'DisplayName', 'Hamming Window');
plot(f, abs(H_hann), 'y', 'DisplayName', 'Hann Window');
plot(f, abs(H_black), 'b', 'DisplayName', 'Blackman Window');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Comparison of Low-Pass Filters Using Different Windows (Non-dB)');
legend;
grid on;
hold off;