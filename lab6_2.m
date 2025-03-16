% Task 2: Band-Pass and Band-Reject Filters
% Original Bandpass Filter Implementation
clc; clear; close all;
fs = 24000; 
fc1 = 1600; 
fc2 = 5000; 
f=(fc1+fc2)/2/fs;
BW = (fc2-fc1)/fs;   
R1 = 1 - 3*BW;
K1 = (1 - 2*R1*cos(2*pi*f) + R1^2) / (2 - 2*cos(2*pi*f));
a0 =1- K1;
a1 = 2*(K1-R1)*cos(2*pi*f);
a2 =R1^2-K1;
b1 = 2*R1*cos(2*pi*f);
b2=-R1^2;
delta=[1,zeros(1,200)];
h(1)= a0* delta(1);
h(2)= a0* delta(2) + a1* delta(1) +b1*h(1);
for i=3:length(delta)
    h(i)=a0* delta(i)+a1*delta(i-1)+a2*delta(i-2)+b1*h(i-1)+b2*h(i-2);
end
H=fft(h);
subplot(321);plot(1:100,h(1:100)); title('Impulse response of 2nd order passband IIR');grid on;
subplot(322);plot(1:100,abs(H(1:100))); title('Frequency response of 2nd order passband IIR');grid on;
%4th order pass-band IIR filter
y(1)=a0*h(1);
y(2)=a0*h(2) + a1*h(1)+ b1* y(1);
for i=3:length(h)
    y(i)=a0* h(i)+a1*h(i-1)+a2*h(i-2)+b1*y(i-1)+b2*y(i-2);
end
Y=fft(y);
subplot(323);plot(1:100,y(1:100)); title('Impulse response of 4th order passband IIR');grid on;
subplot(324);plot(1:100,abs(Y(1:100))); title('Frequency response of 4th order passband IIR');grid on;
t=linspace(0,1,100);
x_noisy =sin(2*pi*100*t) + sin(2*pi*2300*t);
x_filtered = conv(x_noisy,y);
subplot(325);plot(1:length(x_noisy),x_noisy); grid on;title('Noisy signal');
subplot(326);plot(1:100,x_filtered(1:100)); grid on;title('After filterng');

% Band Rejection Filter Implementation (New Figure)
figure;
% Parameters are the same as before
fs = 24000;
fc1 = 1600;
fc2 = 5000;
f = (fc1+fc2)/2/fs;
BW = (fc2-fc1)/fs;

% For band rejection (notch) filter, we need to modify the coefficients
R1 = 1 - 3*BW;
K1 = (1 - 2*R1*cos(2*pi*f) + R1^2) / (2 - 2*cos(2*pi*f));
a0 = K1;                       % Changed for notch filter
a1 = -2*K1*cos(2*pi*f);        % Changed for notch filter
a2 = K1;                       % Changed for notch filter
b1 = 2*R1*cos(2*pi*f);
b2 = -R1^2;

% Original signal for impulse response
delta = [1, zeros(1, 200)];

% 2nd order band rejection filter
h_notch(1) = a0 * delta(1);
h_notch(2) = a0 * delta(2) + a1 * delta(1) + b1 * h_notch(1);
for i = 3:length(delta)
    h_notch(i) = a0 * delta(i) + a1 * delta(i-1) + a2 * delta(i-2) + b1 * h_notch(i-1) + b2 * h_notch(i-2);
end

% Frequency domain analysis
H_notch = fft(h_notch);

% 4th order band rejection filter (cascading two 2nd order filters)
y_notch(1) = a0 * h_notch(1);
y_notch(2) = a0 * h_notch(2) + a1 * h_notch(1) + b1 * y_notch(1);
for i = 3:length(h_notch)
    y_notch(i) = a0 * h_notch(i) + a1 * h_notch(i-1) + a2 * h_notch(i-2) + b1 * y_notch(i-1) + b2 * y_notch(i-2);
end

% Frequency domain analysis for 4th order filter
Y_notch = fft(y_notch);

% Create test signal (same as before)
t = linspace(0, 1, 100);
x_noisy = sin(2*pi*100*t) + sin(2*pi*2300*t);
x_filtered_notch = conv(x_noisy, y_notch);

% Compute frequency responses for notch filtering
X_filtered_notch = fft(x_filtered_notch(1:100));

% Plot 2nd order band rejection filter
subplot(321); plot(1:100, h_notch(1:100)); 
title('Impulse response of 2nd order band rejection IIR'); grid on;
ylabel('Amplitude'); xlabel('Sample');

subplot(322); plot((0:99)/200*fs, abs(H_notch(1:100))); 
title('Frequency response of 2nd order band rejection IIR'); grid on;
ylabel('Magnitude'); xlabel('Frequency (Hz)');

% Plot 4th order band rejection filter
subplot(323); plot(1:100, y_notch(1:100)); 
title('Impulse response of 4th order band rejection IIR'); grid on;
ylabel('Amplitude'); xlabel('Sample');

subplot(324); plot((0:99)/200*fs, abs(Y_notch(1:100))); 
title('Frequency response of 4th order band rejection IIR'); grid on;
ylabel('Magnitude'); xlabel('Frequency (Hz)');


x_filtered_notch = filter([a0 a1 a2], [1 -b1 -b2], x_noisy);

% Plot the time domain signals
subplot(325); 
t = (0:length(x_noisy)-1)/fs; 
plot(t, x_noisy); 
grid on; 
title('Noisy signal (time domain)'); 
ylabel('Amplitude'); 
xlabel('Time (s)');

subplot(326); 
plot(t, x_filtered_notch);  
grid on; 
title('After band rejection filtering (time domain)');
ylabel('Amplitude'); 
xlabel('Time (s)');