n=0:100;
s1=cos(2*pi*0.05*n);
s2=cos(2*pi*0.47*n);
x=s1+s2;
M=input('Enter the order of moving average filter: ');
h=ones(1,M);
y=filter(h, 1,x)/M; %recursive/non-recursive/input signal
clf;
subplot(221);plot(n, s1);grid on;
axis([0, 100, -2, 2]); %first 2 - x, last 2 - y
xlabel('Time(s)'); ylabel('Amplitude'); title('Low frequency signal');

subplot(222);plot(n, s2); grid on;
axis([0, 100, -2, 2]); %first 2 - x, last 2 - y
xlabel('Time(s)'); ylabel('Amplitude'); title('High frequency signal');

subplot(223);plot(n, x); grid on;
axis([0, 100, -2, 2]); %first 2 - x, last 2 - y
xlabel('Time(s)'); ylabel('Amplitude'); title('Input signal');

subplot(224);plot(n, y); grid on;
axis([0, 100, -2, 2]); %first 2 - x, last 2 - y
xlabel('Time(s)'); ylabel('Amplitude'); title('Output signal');
axis;
figure;
X=abs(fft(x));
Y=abs(fft(y));
subplot(211); plot(X(1:50)); title('Magniude of input signal');
subplot(212); plot(Y(1:50)); title('Magniude of output signal');