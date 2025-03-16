clc;clear;close all;
n=6; r=0.5; fc=0.3;
a0=4.18e-3;
a1=-2.51e-02;
a2=6.28e-02;
a3=-8.37e-02;
a4=6.28e-02;
a5=-2.51e-02;
a6=4.18e-03;
b1=-2.31e+00;
b2=-3.29;
b3=-2.9;
b4=-1.69;
b5=-6.02e-01;
b6=-1.03e-01;
a=[a0,a1,a2,a3,a4,a5,a6];
b=[b1, b2, b3, b4, b5, b6];
fvtool(a,b);
t=0:0.001:10;
f1=220; f2=150; f3=420;
s=sin(2*pi*f1*t)+0.7*sin(2*pi*f2*t)+0.3*sin(2*pi*f3*t); %noisy signal
s_out=filter(a,b,s); %after filtering
S=fft(s, 512);
S_out=fft(s_out,512);
subplot(211); plot(t(1:100),s(1:100),'b',t(1:100),s_out(1:100),'r'); grid on; legend('s(t)', 's_out(t)');

subplot(212); plot(f(1:100),abs(S(1:100)),'b',f(1:100),abs(S_out(1:100)),'r'); grid on; legend('S(t)', 'S_out(t)');