% Signal Smoothing by Averaging
clf; clc; clear;
N = 51;
noise = 0.8*(rand(1,N) - 0.5); % Generate random noise-uniform distributed
m = 0:N-1;
s = 2*m.*(0.9.^m); % Generate uncorrupted signal (original)
x_noisy = s + noise; % Generate noise corrupted signal-transpose from col to row
subplot(211);
plot(m,noise,'r', m,s,'b', m,x_noisy,'g');
grid on; title('Before filtering'); xlabel('Time(s)'); ylabel('Amplitude');
legend('Noise', 'Original signal', 'Noisy signal');

x1=[0,0,x_noisy]; x2=[0,x_noisy,0]; x3=[x_noisy,0,0];
y=(x1+x2+x3)/3; %output of filter
subplot(212);
plot(m,y(2:N+1),'r',m,x_noisy,'b');
grid on; title('After filtering'); xlabel('Time(s)'); ylabel('Amplitude');
legend('Filtered signal','Noisy signal');