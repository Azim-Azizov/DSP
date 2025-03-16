%task1
clc;
%creating noisy signal
fs = 1000;
t = 0:1/fs:1-1/fs;
f1 = 50;
f2 = 150;
x = sin(2*pi*f1*t) + sin(2*pi*f2*t);
noise = 0.5*randn(size(t));
x_noisy = x + noise;
figure;
subplot(2,1,1);
plot(t,x);title('Normal Signal'); grid on;xlabel('Time (seconds)');ylabel('Amplitude');

subplot(2,1,2);
plot(t, x_noisy);title('Noisy Signal');xlabel('Time (seconds)');ylabel('Amplitude');


% Low-pass filter
figure;
cutoff_freq = 100; % Low-pass filter cutoff frequency
order = 50;
b = fir1(order, cutoff_freq/(fs/2));
x_filtered = filter(b, 1, x_noisy);

subplot(2,1,1);
plot(t, x_filtered);title('Filtered Signal (Low-Pass)');xlabel('Time (seconds)');ylabel('Amplitude');

[b, a] = butter(4, cutoff_freq/(fs/2));
x_filtered_iir = filter(b, a, x_noisy);

subplot(2,1,2);
plot(t, x_filtered_iir);
title('Filtered Signal (IIR Low-Pass) Butterworth');xlabel('Time (seconds)');ylabel('Amplitude');

% High-pass filter
cutoff_freq_highpass = 100; % High-pass filter cutoff frequency
order = 50;

% FIR High-pass Filter
b_hp = fir1(order, cutoff_freq_highpass/(fs/2), 'high');
x_filtered_hp = filter(b_hp, 1, x_noisy);

% IIR High-pass Filter
[b_hp_iir, a_hp_iir] = butter(4, cutoff_freq_highpass/(fs/2), 'high');
x_filtered_hp_iir = filter(b_hp_iir, a_hp_iir, x_noisy);

% Plotting the high-pass filtered signals
figure;
subplot(2,1,1);
plot(t, x_filtered_hp);title('Filtered Signal (FIR High-Pass)');xlabel('Time (seconds)');ylabel('Amplitude');

subplot(2,1,2);
plot(t, x_filtered_hp_iir);title('Filtered Signal (IIR High-Pass) Butterworth');xlabel('Time (seconds)');ylabel('Amplitude');

% Frequency response of low-pass filter
figure;
freqz(b, 1, 1024, fs);
title('Frequency Response of the Low-Pass FIR Filter');

% Frequency response of high-pass filter
figure;
freqz(b_hp_iir, 1, 1024, fs);
title('Frequency Response of the High-Pass IIR Filter');

% Frequency domain
% Angle phase response abs for magnitude

X_noisy = abs(fft(x_noisy));
X_filtered = abs(fft(x_filtered));
X_original = abs(fft(x));

X_filtered_hp = abs(fft(x_filtered_hp));
X_filtered_hp_iir = abs(fft(x_filtered_hp_iir));

n = length(t);
f = (0:n-1) * fs/n; % Frequency vector

figure;

subplot(2,2,1);
plot(f, X_original);
title('Magnitude Spectrum of Original Signal');xlabel('Frequency (Hz)');ylabel('Magnitude');

subplot(2,2,2);
plot(f, X_noisy);
title('Magnitude Spectrum of Noisy Signal');xlabel('Frequency (Hz)');ylabel('Magnitude');

subplot(2,2,3);
plot(f, X_filtered);
title('Magnitude Spectrum of Filtered Signal');xlabel('Frequency (Hz)');ylabel('Magnitude');

subplot(2,2,4);
plot(f, X_filtered_hp);title('Magnitude Spectrum of High-Pass FIR Filtered Signal');xlabel('Frequency (Hz)');ylabel('Magnitude');
