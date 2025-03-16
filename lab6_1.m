clc;clear; close all;
fc=20;
fs=750;
x=exp(-2*pi*(fc/fs));
%filter coefficients
a_0=1-x;
b_1=x;
%impulse response
delta=[1,zeros(1,1023)];
y(1)=a_0*delta(1);
for i=2:length(delta)
    y(i)=a_0*delta(i)+b_1*y(i-1);
end
%plot(d);title('Impulse response of first order IIR filter');
H1=fft(y);
subplot(321); plot((1:100),y(1:100)); title('Impulse response of first order IIR filter'); grid on;
subplot(322); plot((1:512)/1024,abs(H1(1:512))); title('Frequency response of first order IIR filter'); grid on;
z(1)=a_0*y(1);
for i=2:length(delta)
    z(i)=a_0*y(i)+b_1*z(i-1);
end
H2=fft(z);
subplot(323); plot((1:100),z(1:100)); title('Impulse response of second order IIR filter'); grid on;
subplot(324); plot((1:512)/1024,abs(H2(1:512))); title('Frequency response of second order IIR filter'); grid on;
t=linspace(0,1,100);
x_noisy=sin(2*pi*20*t)+sin(2*pi*300*t);
x(1)=a_0*x_noisy(1); %first sample of noise
for i=2:length(x_noisy)
    x(i)=a_0*x_noisy(i)+b_1*x(i-1); %remove noise
end
subplot(325); plot(1:length(x), x_noisy); title('Noisy signal'); grid on;
subplot(326); plot(1:length(x), x); title('After filtering signal'); grid on;
