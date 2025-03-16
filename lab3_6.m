%Task6
clc; clear; close all;

% Filter parameters
fs = 44100;      
fc1 = 5000;     
fc2 = 15000;     
M = 101;         

% Normalize frequencies
f1 = fc1/fs;
f2 = fc2/fs;

% Time index
n = 1:M;
n_shifted = n - (M+1)/2;  

% Initialize filter coefficient arrays
h_bp = zeros(1,M);
h_bs = zeros(1,M);

% Calculate filter coefficients
for i = 1:M
    if n_shifted(i) == 0 
        h_bp(i) = 2*(f2 - f1);
        h_bs(i) = 1 - 2*(f2 - f1);
    else
        % Band-pass filter
        h_bp(i) = (sin(2*pi*f2*n_shifted(i)) - sin(2*pi*f1*n_shifted(i)))/(pi*n_shifted(i));
        
        % Band-stop filter
        h_bs(i) = (sin(2*pi*f1*n_shifted(i)) - sin(2*pi*f2*n_shifted(i)))/(pi*n_shifted(i));
    end
end

% Apply Hamming window
window = hamming(M)';
h_bp = h_bp .* window;
h_bs = h_bs .* window;

% Calculate frequency responses
[H_bp, f] = freqz(h_bp, 1, 1024, fs);
[H_bs, f] = freqz(h_bs, 1, 1024, fs);

% Plot impulse responses
figure;
subplot(2,1,1);
stem(n, h_bp, 'color', [0.2 0.6 0.8], 'LineWidth', 1.5, 'MarkerFaceColor', [0.2 0.6 0.8]);
title('Band-Pass Filter Impulse Response (Hamming Window)');
xlabel('Samples'); ylabel('Amplitude');
grid on;

subplot(2,1,2);
stem(n, h_bs, 'color', [0.8 0.4 0.2], 'LineWidth', 1.5, 'MarkerFaceColor', [0.8 0.4 0.2]);
title('Band-Stop Filter Impulse Response (Hamming Window)');
xlabel('Samples'); ylabel('Amplitude');
grid on;

% Plot frequency responses
figure;
subplot(2,1,1);
plot(f, abs(H_bp), 'color', [0.4 0.8 0.2], 'LineWidth', 2); hold on;
plot(f, abs(H_bs), 'color', [0.9 0.3 0.3], 'LineWidth', 2);
title('Frequency Response (Linear Scale)');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
legend('Band-Pass Filter', 'Band-Stop Filter');
grid on; xlim([0 fs/2]);
hold off;

subplot(2,1,2);
plot(f, 20*log10(abs(H_bp)), 'color', [0.4 0.8 0.2], 'LineWidth', 2); hold on;
plot(f, 20*log10(abs(H_bs)), 'color', [0.9 0.3 0.3], 'LineWidth', 2);
title('Frequency Response (dB Scale)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend('Band-Pass Filter', 'Band-Stop Filter');
grid on; xlim([0 fs/2]); ylim([-120 10]);
hold off;