clc; clear; close all;

Wp = [0.2 0.4];

[b, a] = cheby1(2, 0.3, Wp, "bandpass");

%Impulse response
[h, t]=impz(b, a, 1024);

%Frequency response
[G, f]=freqz(b, a, 512);    

figure;
subplot(311); plot(t, h); title("Impulse response of the filter"); grid on;
subplot(312); plot(f, abs(G)); title("Frequency response of the filter (non-dB)"); grid on;
subplot(313); plot((0 : 511) / 1024, abs(10*log10(G(1 : 512)))); title("Frequency response of the filter (dB)"); grid on;

%Applying filter to an input signal
fs = 1000;
t2 = 0:1/fs:1;
f1 = 15; f2 =460; f3 = 220; f4 = 250;
x = sin(2*pi*f1*t2) + sin(2*pi*f2*t2) + sin(2*pi*f3*t2) + 2*sin(2*pi*f4*t2);

y  = filter(b, a, x);

figure;
subplot(221); plot(t2, x); title("Input signal"); grid on;
subplot(222); plot(t2,y); title("Filtered input signal"); grid on;

X_freq = fft(x, 512);
Y_freq = fft(y, 512);

subplot(223); plot((1: 512), abs(X_freq(1 : 512))); title("Spectrum before filtering"); grid on;
subplot(224); plot((1: 512), abs(Y_freq(1 : 512))); title("Spectrum after filtering"); grid on;