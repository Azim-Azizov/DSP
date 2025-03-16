fs = 1000; 
t = 0:1/fs:1; % Time vector 
f = 5; % Signal frequency 
x = sin(2*pi*f*t); 

figure;
plot(t, x);
xlabel('Time (s)');
ylabel('Amplitude');
title('Sinusoidal Signal');
grid on;

f = 5; % Frequency of the square wave
sq_wave = square(2 * pi * f * t);

figure;
plot(t, sq_wave);
xlabel('Time (s)');
ylabel('Amplitude');
title('Square Wave Signal');
grid on;

%2 Samplig and Aliasing
f_signal = 10; 
fs_low = 15; % Low sampling rate 
fs_high = 50; % High sampling rate

t_low = 0:1/fs_low:1;
t_high = 0:1/fs_high:1;

x_low = sin(2*pi*f_signal*t_low);
x_high = sin(2*pi*f_signal*t_high);

figure;
subplot(2,1,1);
stem(t_low, x_low, 'r'); hold on;
plot(t_high, x_high, 'b');
xlabel('Time (s)'); ylabel('Amplitude');
title('Aliasing Effect (Red: Low Sampling, Blue: Original Signal)');
legend('Low Sampling', 'High Sampling');
grid on;

%on square wave 

f_signal = 10; % Signal frequency
fs_low = 15; % Low sampling rate (below Nyquist rate)
fs_high = 150; % High sampling rate

t_low = 0:1/fs_low:1;
t_high = 0:1/fs_high:1;

x_low = square(2*pi*f_signal*t_low);
x_high = square(2*pi*f_signal*t_high);

subplot(2,1,2);
stem(t_low, x_low, 'r'); hold on;
plot(t_high, x_high, 'b');
xlabel('Time (s)'); ylabel('Amplitude');
title('Aliasing Effect (Red: Low Sampling, Blue: Original Signal)');
legend('Low Sampling', 'High Sampling');
grid on;

%fft 
N = length(x);
X = fft(x);
f_axis = linspace(-fs/2, fs/2, N);

figure;
plot(f_axis, abs(fftshift(X)));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of the Signal');
grid on;





