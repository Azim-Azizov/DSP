

























































































































































































































% Lab 8
%% Task 1 - Lowpass FIR/IIR comparison
fs=8000;
t=0:1/fs:1;
signal=sin(2*pi*500*t)+sin(2*pi*2000*t);
noise=0.5*randn(size(t));
signal=signal+noise;
filtered_FIR=filter(Num,1,signal);
filtered_IIR=filtfilt(SOS,G,signal);
figure;
subplot(311);plot(t(1:100),signal(1:100));grid on;
title('Original signal');xlabel('Time (s)');ylabel('Amplitude');

subplot(312);plot(t(1:100),filtered_FIR(1:100));grid on;
title('Filtered signal (FIR)');xlabel('Time (s)');ylabel('Amplitude');

subplot(313);plot(t(1:100),filtered_IIR(1:100));grid on;
title('Filtered signal (IIR)');xlabel('Time (s)');ylabel('Amplitude');

sgtitle('Comparison of FIR and IIR filtering');

%Compute FFT
N=length(signal);
f=linspace(0,fs/2,N/2);
X_original=abs(fft(signal));
X_FIR=abs(fft(filtered_FIR));
X_IIR=abs(fft(filtered_IIR));

%Normalize and keep only positive frequencies
X_original=X_original(1:N/2)/max(X_original);
X_FIR=X_FIR(1:N/2)/max(X_FIR);
X_IIR=X_IIR(1:N/2)/max(X_IIR);

figure;
subplot(311);plot(f,X_original);grid on;
title('Frequency response of original signal');xlabel('Frequency (Hz)');ylabel('Magnitude');

subplot(312);plot(f,X_FIR);grid on;
title('Frequency response of filtered signal (FIR)');xlabel('Frequency (Hz)');ylabel('Magnitude');

subplot(313);plot(f,X_IIR);grid on;
title('Frequency response of filtered signal (IIR)');xlabel('Frequency (Hz)');ylabel('Magnitude');

sgtitle('Comparison of frequency responses');

%% Task 2 - Highpass FIR/IIR comparison
fs=8000;
t=0:1/fs:1;
signal=sin(2*pi*500*t)+sin(2*pi*2500*t);
noise=0.5*randn(size(t));
signal=signal+noise;
filtered_FIR=filter(Num1,1,signal);
filtered_IIR=filtfilt(SOS1,G1,signal);
figure;
subplot(311);plot(t(1:100),signal(1:100));grid on;
title('Original signal');xlabel('Time (s)');ylabel('Amplitude');

subplot(312);plot(t(1:100),filtered_FIR(1:100));grid on;
title('Filtered signal (FIR)');xlabel('Time (s)');ylabel('Amplitude');

subplot(313);plot(t(1:100),filtered_IIR(1:100));grid on;
title('Filtered signal (IIR)');xlabel('Time (s)');ylabel('Amplitude');

sgtitle('Comparison of FIR and IIR filtering');

%Compute FFT
N=length(signal);
f=linspace(0,fs/2,N/2);
X_original=abs(fft(signal));
X_FIR=abs(fft(filtered_FIR));
X_IIR=abs(fft(filtered_IIR));

%Normalize and keep only positive frequencies
X_original=X_original(1:N/2)/max(X_original);
X_FIR=X_FIR(1:N/2)/max(X_FIR);
X_IIR=X_IIR(1:N/2)/max(X_IIR);

figure;
subplot(311);plot(f,X_original);grid on;
title('Frequency response of original signal');xlabel('Frequency (Hz)');ylabel('Magnitude');

subplot(312);plot(f,X_FIR);grid on;
title('Frequency response of filtered signal (FIR)');xlabel('Frequency (Hz)');ylabel('Magnitude');

subplot(313);plot(f,X_IIR);grid on;
title('Frequency response of filtered signal (IIR)');xlabel('Frequency (Hz)');ylabel('Magnitude');

sgtitle('Comparison of frequency responses');

%% Task 3 - Bandpass FIR/IIR comparison
fs=8000;
t=0:1/fs:1;
signal = sin(2*pi*100*t) + sin(2*pi*800*t) + sin(2*pi*370*t);
noise=0.5*randn(size(t));
signal=signal+noise;
filtered_FIR=filter(Num2,1,signal);
filtered_IIR=filtfilt(SOS2,G2,signal);
figure;
subplot(311);plot(t(1:100),signal(1:100));grid on;
title('Original signal');xlabel('Time (s)');ylabel('Amplitude');

subplot(312);plot(t(1:100),filtered_FIR(1:100));grid on;
title('Filtered signal (FIR)');xlabel('Time (s)');ylabel('Amplitude');

subplot(313);plot(t(1:100),filtered_IIR(1:100));grid on;
title('Filtered signal (IIR)');xlabel('Time (s)');ylabel('Amplitude');

sgtitle('Comparison of FIR and IIR filtering');

%Compute FFT
N=length(signal);
f=linspace(0,fs/2,N/2);
X_original=abs(fft(signal));
X_FIR=abs(fft(filtered_FIR));
X_IIR=abs(fft(filtered_IIR));

%Normalize and keep only positive frequencies
X_original=X_original(1:N/2)/max(X_original);
X_FIR=X_FIR(1:N/2)/max(X_FIR);
X_IIR=X_IIR(1:N/2)/max(X_IIR);

figure;
subplot(311);plot(f,X_original);grid on;
title('Frequency response of original signal');xlabel('Frequency (Hz)');ylabel('Magnitude');

subplot(312);plot(f,X_FIR);grid on;
title('Frequency response of filtered signal (FIR)');xlabel('Frequency (Hz)');ylabel('Magnitude');

subplot(313);plot(f,X_IIR);grid on;
title('Frequency response of filtered signal (IIR)');xlabel('Frequency (Hz)');ylabel('Magnitude');

sgtitle('Comparison of frequency responses');


% Lab 9
%% Task 1 
a=audioinfo('omni_drum.mp3');
disp(a);
whos a;

a.Title='My Song';
a.Rating=[10, 5];
mystructure=a;
mystructure.b=5;

sampleRate=a.TotalSamples/a.Duration;
disp(['Calculated Sample Rate:', num2str(sampleRate)]);

%% Task 2 - creating sinewave
 %Create a sinewave with a frequency of 1000 Hz, amplitude of .2 and phase of  15°. 
% What would be the minimum sampling rate that this sinewave would require? 
fs=3000;
t=0:1/fs:2;
f=1000;
A=.5;
w=15*pi/180;
y=A*sin(2*pi*f*t+w);
sound(y,fs,16);
plot(t,y);xlabel('Time (s)');ylabel('Amplitude');
title('Sinewave with a frequency of 1000 Hz, amplitude of 0.2 and phase of 15 degree');


%% Task 3 - Phase cancellation
fs=44100;
f1=440;
f2=440;
A1=.3;
A2=.3;
t=0:1/fs:5;
w1=0*pi/180;
w2=180*pi/180;
y1=A1*sin(2*pi*f1*t+w1);
y2=A2*sin(2*pi*f2*t+w2);
y=(y1+y2)/2;
sound(y,fs,16);
plot(t,y);title('Phase cancellation');xlabel('Time (s)');ylabel('Amplitude');

% Change one frequency to 880 Hz
f2=880;
y2=A2*sin(2*pi*f2*t+w2);
y=(y1+y2)/2;
sound(y,fs,16);
figure;
subplot(221);plot(t,y);
title('Changed frequency to 880 Hz');xlabel('Time (s)');ylabel('Amplitude');

% Change one frequency to 441 Hz
f2=441;
y2=A2*sin(2*pi*f2*t+w2);
y=(y1+y2)/2;
sound(y,fs,16);
subplot(222);plot(t,y);
title('Changed frequency to 441 Hz');xlabel('Time (s)');ylabel('Amplitude');

% Phase changed to 179 degrees
w2=179*pi/180;
y2=A2*sin(2*pi*f2*t+w2);
sound(y,fs,16);
subplot(223);plot(t,y);
title('Phase changed to 179 degrees');xlabel('Time(s)');ylabel('Amplitude');

% Phase changed to 181 degrees
w2=181*pi/180;
y2=A2*sin(2*pi*f2*t+w2);
sound(y,fs,16);
subplot(224);plot(t,y);
title('Phase changed to 181 degrees');xlabel('Time(s)');ylabel('Amplitude');

%% Task 4 - Binaural beats/Cocktail party effect
fs=44100;
t=0:1/fs:5;
f1=300;
f2=310;
w=0*pi/180;
A=.5;
y1=A*sin(2*pi*f1*t+w);
y2=A*sin(2*pi*f2*t+w);
y=[y1;y2];
sound(y,fs,16);

% Change frequency difference to 2 Hz
f2=f1+2;
y2=A*sin(2*pi*f2*t+w);
y=[y1;y2];
sound(y,fs,16);

% Play 1000 Hz and 1010 Hz
f1=1000;
f2=1010;
y1=A*sin(2*pi*f1*t+w);
y2=A*sin(2*pi*f2*t+w);
y=[y1;y2];
sound(y,fs,16);

% Create cocktail party effect
f1=300;
f2=330;
y1=A*sin(2*pi*f1*t+w);
y2=A*sin(2*pi*f2*t+w);
y=[y1;y2];
sound(y,fs,16);

%% Task 5 - Complex tones via additive synthesis
fs=44100;
t=0:1/fs:5;
f1=400; f2=2*f1; f3=3*f1; f4=4*f1;
A1=.3; A2=A1/2; A3=A1/3; A4=A1/4;
w=0;
y1=A1*sin(2*pi*f1*t+w);
y2=A2*sin(2*pi*f2*t+w);
y3=A3*sin(2*pi*f3*t+w);
y4=A4*sin(2*pi*f4*t+w);
y=(y1+y2+y3+y4)/4;
sound(y,fs,16);

%% Task 6 - Melodies
fs=44100;
t=0:1/fs:5;
A1=.3; A2=A1/2; A3=A1/3; A4=A1/4;
f1=440; f2=f1*2; f3=f1*3; f4=f1*4;
w=0;
y1=A1*sin(2*pi*f1*t+w);
y2=A2*sin(2*pi*f2*t+w);
y3=A3*sin(2*pi*f3*t+w);
y4=A4*sin(2*pi*f4*t+w);
y=[y1 y2 y3 y4];
soundsc(y,fs);

x=y(1:2:end);
soundsc(x,fs);

%% Task 7 - Modulation
fs=44100;
t=0:1/fs:5;
f=440;
y=sin(2*pi*f*t);
% frequency modulation
fm=modulate(y,20,fs,'fm');
soundsc(fm,fs);
% amplitude modulation
am=modulate(y,.5,fs,'am');
soundsc(am,fs);
% phase modulation
pm=modulate(y,20,fs,'pm');
soundsc(pm,fs);



%% Task 8 - FFT
%Create a signal y at 440 Hz and plot its spectrum. Use an adequate axis scaling 
%(help axis). 
% Plot the spectrum of a complex sinewave at different amplitudes. 
% Do the plots include phase information of the signal? Explain.
fs=44100;
t=0:1/fs:5;
f=440;
y=sin(2*pi*f*t);
N=fs;
Y=fft(y,N)/N;
magTransform=abs(Y);
faxis=linspace(-fs/2,fs/2,N);
plot(faxis, fftshift(magTransform));
xlabel('Frequency (Hz)');
axis([0 2000 0 max(magTransform)]);

%% Task 9 - Spectrogram
%Create a quadratic chirp starting at 0 Hz and reaching 440 Hz at 1 s using c = 
%chirp(t,0,1,440,'q'). Plot its spectrogram. 
% Plot a spectrogram of the melody that you created before. 
% Try the effect of using different window lengths on the spectrogram. What 
%happens when you modify this length?
fs=44100;
t=0:1/fs:5;
f=440;
y=sin(2*pi*f*t);
win=128;
hop=win/2;
nfft=win;
spectrogram(y,win,hop,nfft,fs,'yaxis');

% quadratic chirp
c=chirp(t,0,1,440,'q');
spectrogram(c,win,hop,nfft,fs,'yaxis');

%% Task 10 - Shephard Tones
fs=16000;
d=1;
t=0:1/fs:d-1/fs;
fmin=300;
fmax=3400;
n=12;
l=mod(([0:n-1]/n)'*ones(1,fs*d)+ones(n,1)*(t/(d*n)),1);
f=fmin*(fmax/fmin).^l;
p=2*pi*cumsum(f,2)/fs;
p=diag((2*pi*floor(p(:,end)/(2*pi)))./p(:,end))*p;
s=sin(p);
a=0.5-0.5*cos(2*pi*l);
w=sum(s.*a)/n;
w=repmat(w,1,3);
specgram(w,2048,fs,2048,1800);
ylim([0 4000]);
soundsc(w,fs);
audiowrite('shephard.wav',w,fs);

%% Task 11 - Amplitude modulation and Ring modulation
load handel;
index=1:length(y);
Fc=5;
A=.5;
w=0;
trem=(w*pi/180+A*sin(2*pi*index*(Fc/Fs)))';
y=y.*trem;
soundsc(y,Fs);

Fc=1000;
trem=(w*pi/180+A*sin(2*pi*index*(Fc/Fs)))';
y=y.*trem;
soundsc(y,Fs);

%% Task 12 - Filter design
%Listen to the result. Listen to the original tone and compare both. 
% Plot the beginning of both versions of the tone 
%with plot(o(1:800)),figure,plot(y(1:800)). Can you find any differences? 
% Try modifying the filter cutoff and the filter order. How do they affect the sound? 
% Repeat the procedure but using a highpass filter (help designfilt)
fs=44100;
t=0:1/fs:5;
y=sin(2*pi*500*t)+sin(2*pi*2500*t);
noise=0.5*randn(size(t));
y=y+noise;
% Lowpass filter cutoff 442 order 10
cutoff=442/(fs/2);
order=10;
d=designfilt('lowpassfir','CutoffFrequency',cutoff,'FilterOrder',order);
o=filter(d,y);
soundsc(o,fs);

figure;
subplot(221);plot(y(1:800));xlabel('Sample Index');ylabel('Amplitude');
title('Original signal (first 800 samples)');

subplot(222);plot(o(1:800));xlabel('Sample Index');ylabel('Amplitude');
title('Filtered signal (Lowpass, cutoff 442 Hz, order 10)');

% Change cutoff and filter order
cutoff=600/(fs/2);
order=20;
d=designfilt('lowpassfir','CutoffFrequency',cutoff,'FilterOrder',order);
o2=filter(d,y);
soundsc(o2,fs);

subplot(223);plot(o2(1:800));xlabel('Sample Index');ylabel('Amplitude');
title('Filtered signal (Lowpass, cutoff 600 Hz, order 20)');

% Highpass filter cutoff 420 order 10
cutoff=442/(fs/2);
order=10;
d_high=designfilt('highpassfir','CutoffFrequency',cutoff,'FilterOrder',order);
o_high=filter(d_high,y);
soundsc(o_high,fs);

subplot(224);plot(o_high(1:800));xlabel('Sample Index');ylabel('Amplitude');
title('Filtered signal (Highpass, cutoff 442 Hz, order 10)');

% Lab 10
%% Task 1
% cceps function to show an echo in the signal and remove echo from signal. - 
%Generate a sine of frequency 45 Hz, sampled at 100 Hz. Add an echo with half the 
%amplitude and 0.2 s later.  Play (sound) given signal. 
%- Compute the complex cepstrum of the signal. Notice the echo at 0.2 s. 
%- Remove echo from signal using necessary filter type. 
%- Plot obtained signal signal in time anf frequency domain 
%- Play obtained signal. 
%- Load into program another audio signal. Add an echo with 0.4 amplitude and 0.3 s 
%delay.  
%Digital Signal and Data Processing - Eliminate echo from this signal using homomorphic analysis. 
f=45;
fs=100;
duration=2;
t=0:1/fs:duration;
x=sin(2*pi*f*t);
echo_amplitude=.5;
echo_delay=.2;
delay_samples=round(echo_delay*fs);
echo=[zeros(1,delay_samples), echo_amplitude*x];
x_original=[x, zeros(1, delay_samples)];
y=x_original+echo;

sound(y,8000);
pause(duration+echo_delay+0.5);

figure;
subplot(221);plot((0:length(x_original)-1)/fs,x_original);
title('Original signal');xlabel('Time (s)');ylabel('Amplitude');

X=abs(fft(x_original));
faxis=linspace(0,fs,length(X));
subplot(222);plot(faxis,X);
title('Original signal in frequency domain');
xlabel('Frequency (Hz)');ylabel('Magnitude');

subplot(223);plot((0:length(y)-1)/fs,y);
title('Signal with echo');xlabel('Time (s)');ylabel('Amplitude');

Y=abs(fft(y));
faxis=linspace(0,fs,length(Y));
subplot(224);plot(faxis,Y);
title('Signal with echo in frequency domain');
xlabel('Frequency (Hz)');ylabel('Magnitude');

% compute complex cepstrum
[xhat, nd]=cceps(y);

figure;
plot((0:length(xhat)-1)/fs,xhat);
title('Signal in complex cepstrum');xlabel('Quefrency (s)');ylabel('Amplitude');

% remove echo
cepstrum_cleaned=xhat;
cepstrum_cleaned(delay_samples)=0;
y_clean=icceps(cepstrum_cleaned, nd);

sound(real(y_clean),8000);
pause(duration+echo_delay+0.5);

figure;
subplot(211);plot((0:length(y_clean)-1)/fs,real(y_clean));
title('Signal after echo removal');xlabel('Time (s)');ylabel('Amplitude');

Y_clean=abs(fft(real(y_clean)));
faxis=linspace(0,fs,length(Y_clean));
subplot(212);plot(faxis,Y_clean);
title('Signal after echo removal in frequency domain');
xlabel('Frequency (Hz)');ylabel('Magnitude');

%% Task 1.2
fs=44100;
duration=5;
t=0:1/fs:duration;
x=sin(2*pi*500*t)+sin(2*pi*2500*t);
noise=0.5*randn(size(t));
x=x+noise;
echo_amplitude=.5;
echo_delay=.2;
delay_samples=round(echo_delay*fs);
echo=[zeros(1,delay_samples), echo_amplitude*x];
x_original=[x, zeros(1, delay_samples)];
y=x_original+echo;

sound(y,8000);
pause(duration+echo_delay+0.5);

figure;
subplot(221);plot((0:length(x_original)-1)/fs,x_original);
title('Original signal');xlabel('Time (s)');ylabel('Amplitude');

X=abs(fft(x_original));
faxis=linspace(0,fs,length(X));
subplot(222);plot(faxis,X);
title('Original signal in frequency domain');
xlabel('Frequency (Hz)');ylabel('Magnitude');

subplot(223);plot((0:length(y)-1)/fs,y);
title('Signal with echo');xlabel('Time (s)');ylabel('Amplitude');

Y=abs(fft(y));
faxis=linspace(0,fs,length(Y));
subplot(224);plot(faxis,Y);
title('Signal with echo in frequency domain');
xlabel('Frequency (Hz)');ylabel('Magnitude');

% compute complex cepstrum
[xhat, nd]=cceps(y);

figure;
plot((0:length(xhat)-1)/fs,xhat);
title('Signal in complex cepstrum');xlabel('Quefrency (s)');ylabel('Amplitude');

% remove echo
cepstrum_cleaned=xhat;
cepstrum_cleaned(delay_samples)=0;
y_clean=icceps(cepstrum_cleaned, nd);

sound(real(y_clean),8000);
pause(duration+echo_delay+0.5);

figure;
subplot(211);plot((0:length(y_clean)-1)/fs,real(y_clean));
title('Signal after echo removal');xlabel('Time (s)');ylabel('Amplitude');

Y_clean=abs(fft(real(y_clean)));
faxis=linspace(0,fs,length(Y_clean));
subplot(212);plot(faxis,Y_clean);
title('Signal after echo removal in frequency domain');
xlabel('Frequency (Hz)');ylabel('Mgnitude');

Lec 7:
Analog vs. Digital Filters
•	Digital signals often originate from analog electronics, raising the question of whether to filter before or after digitization.
•	While analog filters can be constructed with op amps, resistors, and capacitors, digital filters digitize the analog signal and process it digitally.
Analog vs. Digital Filters (frequency response)
•	Digital filters generally offer better performance in passband ripple, roll-off, and stopband attenuation compared to analog filters.
•	Analog filters have limitations in flatness due to the accuracy of resistors and capacitors, while digital filters are limited by round-off error, making them much flatter.
•	Digital filters exhibit superior roll-off and stopband attenuation, and even with improvements, analog filters cannot match this performance.
Analog vs. Digital Filters (step response)
•	Digital filters demonstrate better step response symmetry (linear phase) compared to analog filters, which have a non-linear phase.
•	Analog filters may exhibit overshoot on one side of the step response, while digital filters overshoot on both sides. Since both are bad
Advantages of Analog Filters
•	Analog filters are still relevant in applications where speed is critical.
•	Analog circuits offer a wider dynamic range in both amplitude and frequency compared to digital systems.
•	Analog filters can handle a broader range of frequencies simultaneously, while digital systems may face overflow issues with large datasets.
Windowed-Sinc vs. Chebyshev
•	Both windowed-sinc (FIR) and Chebyshev (IIR) filters are designed for frequency band separation.
•	A fair comparison requires considering that the Chebyshev's frequency response changes with the cutoff frequency.
•	Windowed-sinc filters generally have better stopband attenuation than Chebyshev filters.
Windowed-Sinc vs. Chebyshev (perfomance)
•	For maximum performance, windowed-sinc filters outperform Chebyshev filters, particularly in scenarios requiring isolation of signals with close frequencies.
•	There are limits to the maximum performance that recursive filters can provide, whereas the windowed-sinc can be pushed to incredible levels.
Windowed-Sinc vs. Chebyshev (speed)
•	IIR filters (e.g., Chebyshev) are generally faster than FIR filters (e.g., windowed-sinc) of comparable performance.
•	The windowed-sinc execution time rises at low and high frequencies because the filter kernel must be made longer to keep up with the greater performance of the recursive filter at these frequencies.
Moving Average vs. Single Pole
•	This section compares time domain filters, specifically a nine-point moving average filter and a single pole recursive filter.
•	Both filters have poor frequency responses, indicating that frequency separation is not their primary use.
•	The moving average filter provides a more rapid step response, while the recursive filter's response is smoother.
Moving Average vs. Single Pole (perfomance)
•	The choice between moving average and single pole filters often depends on the trade-off between development time and execution time.
•	When reducing development time is a priority, the single pole recursive filter is preferred due to its ease of programming and faster execution.
•	When speed is critical, the moving average filter is the better choice, particularly when implemented with recursion and using integers for maximum speed.


Lec 8:
Direction and Distance to the Sound Source
•	Human hearing can determine sound direction but not easily the distance.
•	High frequencies indicate nearby sounds, while low frequencies indicate distant sounds.
•	Echo intervals can provide clues about the size of the listening environment.
Active Sonar
•	Some species use active sonar to locate objects.
•	Bats and dolphins use clicks and squeaks.
•	Some humans can also use active echo localization.
Sound Wave Parts
•	Continuous sounds are perceived through loudness, pitch, and timbre.
•	Loudness measures sound wave intensity, and pitch is the fundamental frequency.
Timbre
•	Timbre is determined by the harmonic content of the signal.
•	Hearing is based on the amplitude of frequencies and is insensitive to their phase.
•	Sound propagation through complex environments randomizes the phase of audio signals.
•	Timbre results from the ear detecting harmonics. 
•	A waveform has only one timbre, but a timbre can have infinite waveforms.
•	Combinations of sine waves with harmonic frequencies sound natural and pleasant.





Lec 9:
Sound Quality vs. Data Rate
•	Digital audio system design balances sound quality and tolerable data rate.
•	The categories include high fidelity music, telephone communication, and compressed speech.
Categories of Digital Audio System
•	High fidelity systems prioritize sound quality, accepting almost any data rate.
•	Telephone communication systems aim for natural sounding speech with a low data rate. (for reducing the system cost)
•	Compressed speech systems prioritize reducing the data rate, tolerating some unnaturalness in sound quality. (Military communication, cellular telephones, and digitally stored speech for voice mail and multimedia)
Sound Quality vs. Data Rate
•	High fidelity music systems use a sampling rate of 44.1 kHz with 16 bits, resulting in a data rate of 706k bits/sec.
•	Natural sounding speech requires about 3.2 kHz bandwidth, whereas music requires 20 kHz.
•	Telecommunication systems typically use an 8 kHz sampling rate for natural sounding speech.
•	FM radio stations broadcast with a bandwidth of almost 20 kHz, while AM radio stations are limited to about 3.2 kHz.
•	Voice-only systems reduce precision to 12 bits or even 8 bits per sample using companding.
•	An 8 kHz sampling rate with 8 bits per sample results in a 64k bits/sec data rate for natural sounding speech.
•	Speech requires less than 10% of the data rate of high fidelity music.
•	Techniques for lowering the data rate are based on compressing the data stream by removing redundancies in speech signals.
•	Linear Predictive Coding (LPC) can reduce the data rate to as little as 2-6k bits/sec, depending on required speech quality.

 
High Fidelity Audio
•	High fidelity audio systems are designed to exceed the limits of human hearing to ensure original music reproduction.
•	Digital audio was popularized by the compact laser disc (CD), which offered superior sound quality compared to older systems.
•	CD stores digital information as dark pits burned on the surface with a laser.
•	An optical sensor detects reflective or nonreflective surfaces during playback, generating binary information.
•	The CD stores about 1 bit per 𝜇𝑚^2, corresponding to 1 million bits per 𝑚𝑚^2, and 15 billion bits per disk.
•	Light cannot be focused to smaller than about one-half wavelength, or 0.3 μm.
•	The raw data rate of a CD playback system is 4.3 million bits per second.
•	Each pit must be no shorter than 0.8 μm, and no longer than 3.5 μm to reduce the error rate due to the optical pickup.
•	Eight-to-fourteen modulation (EFM) is used to ensure binary data complies with bunching requirements.
•	8 bits are converted to 14 bits using a look-up table for storage on the disc, then converted back during playback.
•	Data is encoded using two-level Reed-Solomon coding, combining stereo channels with error detection and correction data.
•	Digital errors are either corrected using redundant data, interpolated between adjacent samples, or muted.
•	These encoding schemes result in the data rate being tripled, from 1.4 Mbits/sec to 4.3 Mbits/sec.
•	After decoding, audio signals are represented as 16 bit samples at 44.1 kHz.
•	Multirate techniques, like converting the digital data to a higher sampling rate (e.g., 176.4 kHz), are commonly used before DAC.
•	Interpolation involves adding three samples with a value of zero between the original samples.
•	Increasing the sample rate results in a smoother signal generated by the DAC.
•	The analog filter only needs to pass frequencies below 20 kHz and block frequencies above 88.2 kHz.
•	Since there are four times as many samples, the number of bits per sample can be reduced from 16 bits to 14 bits, without degrading the sound quality.
•	The sin(𝑥)/𝑥 correction can be part of either the analog or digital filter.

High Fidelity Audio - Stereo
•	Audio systems with more than one channel are said to be in stereo.
•	Multiple channels send sound to the listener from different directions, providing a more accurate reproduction of the original music.
•	Good stereo reproduction makes the listener feel as if the musicians are only a few feet away.
•	High fidelity music has used two channels (left and right) since the 1960s, while motion pictures have used four channels (left, right, center, and surround).
•	Mix-down is an art, aimed at providing the listener with the perception of being there.
•	Four-channel sound in motion pictures is called Dolby Stereo, with the home version called Dolby Surround Pro Logic.
•	The four channels are encoded into the standard left and right channels.
•	The center channel reproduces speech and other visually connected sounds.
•	The surround speakers are placed to the left and right of the listener.
•	The surround channel only contains midrange frequencies and is delayed by 15 to 30 milliseconds.
•	The listener's mind interprets the delayed signal as a reflection from the walls, and ignores it.

Companding
•	Companding reduces the data rate of audio signals by making quantization levels unequal.
•	The loudest sound that can be tolerated (120 dB SPL) is about one million times the amplitude of the weakest sound that can be detected (0 dB SPL).
•	The ear can only distinguish about 120 different loudness levels.
•	If quantization levels are equally spaced, 12 bits are needed for telephone quality speech.
•	Only 8 bits are required if the quantization levels are made unequal, matching human hearing.
•	Small signals require closely spaced levels, while larger signals can use larger spacing.
•	Companding can be carried out by running the analog signal through a nonlinear circuit, using an ADC with unequally spaced steps, or using a linear ADC followed by a digital look-up table.
•	Each option requires the same nonlinearity in a different place.
•	The μ-law algorithm is used in 8-bit digital telecommunication systems in North America and Japan, while A-law is used in Europe.
•	Both use a logarithmic nonlinearity.
•	The curves for μ-law and A-law are nearly identical, with the only significant difference near the origin.
•	μ-law is a smooth curve, and "A" law switches to a straight line.

Speech Synthesis and Recognition
•	Computer generation and recognition of speech are difficult problems.
•	Most commercial products play back digitally recorded segments from a human speaker.

Model of human speech production
•	Nearly all techniques for speech synthesis and recognition are based on the model of human speech production.
•	Most human speech sounds can be classified as either voiced or fricative.
•	Voiced sounds occur when air is forced from the lungs, through the vocal cords, and out of the mouth and/or nose.
•	Voiced sounds are represented by a pulse train generator, with adjustable pitch.
•	Fricative sounds originate as random noise, not from vibration of the vocal cords.
•	This occurs when the air flow is nearly blocked by the tongue, lips, and/or teeth.
•	Fricatives are represented by a noise generator.
•	Both voiced and fricative sound sources are modified by acoustic cavities.
•	Sound propagation through these structures can be represented as a linear filter.
•	Resonance peaks are called the format frequencies.

Spectrogram
•	The voice spectrogram displays speech signals by breaking the audio into short segments and using the FFT to find the frequency spectrum of each segment.
•	These spectra are placed side-by-side and converted into a grayscale image.
•	The segment length is a tradeoff between frequency resolution and time resolution.
•	Voiced sounds have a periodic time domain waveform and a frequency spectrum that is a series of regularly spaced harmonics.
•	Fricatives have a noisy time domain signal and a noisy spectrum.

Speech Synthesis
•	Over a short period, a speech signal can be approximated by specifying excitation (periodic or random noise), frequency of the periodic wave, and coefficients of the digital filter.
•	Continuous speech can be synthesized by continually updating these three parameters about 40 times a second.
•	This approach was used in the Speak & Spell, an electronic learning aid for children.

•	This is also the basis for the linear predictive coding (LPC) method of speech compression.
•	Digitally recorded human speech is broken into short segments, and each is characterized according to the three parameters of the model.
•	The segment information is transmitted or stored as needed, and then reconstructed with the speech synthesizer.

Speech Recognition
•	Speech recognition algorithms try to recognize patterns in extracted parameters by comparing segment information with templates of previously stored sounds.
•	This method does not work very well and is far below the capabilities of human listeners.
•	Words are recognized by their sounds, but also by the context of the sentence, and the expectations of the listener.
•	Listeners hear the correct words for the context, even if exactly the same sounds were produced.
•	This is because of accumulated knowledge about the world.
•	Most speech recognition algorithms rely only on the sound of the individual words, and not on their context.
•	They attempt to recognize words, but not to understand speech.

Disadvantages of Speech Recognition
•	Recognized speech must have distinct pauses between the words.
•	The vocabulary is often limited to only a few hundred words.
•	The algorithm must be trained on each speaker.
Speech Recognition
•	Speech is the quickest and most efficient way for humans to communicate.
•	Speech recognition has the potential of replacing writing, typing, keyboard entry, and electronic control.
•	Progress in speech recognition is related to artificial intelligence, neural networks, and DSP.

Lec 10:
Nonlinear Audio Processing
•	Digital filtering improves audio signals, using techniques like Wiener filtering and deconvolution.
•	Nonlinear techniques are useful for audio processing.
Reducing Wideband Noise in Speech Signals
•	Time varying Wiener Filter reduces wideband noise like tape hiss and wind noise.
•	Linear filtering is ineffective due to frequency overlap between noise and voice signals (200 Hz to 3.2 kHz).
•	Separation is achieved by analyzing the amplitude of each frequency component.
Overlap-add Method
•	Signals and noise can be separated by amplitude; large amplitudes are retained as signal, small amplitudes are discarded as noise. Mid-size frequency components are adjusted in some smooth manner between the two extremes.
•	This technique functions as a time-varying Wiener filter, recalculating frequency response for each segment based on its spectrum. In other words, the filter's frequency response changes from segment-to-segment, as determined by the characteristics of the signal itself.
•	Nonlinear techniques complicate filtering long signals since the overlap-add method becomes invalid due to changing frequency responses, which misalign time-domain segments. A common solution is to split the signal into overlapping segments, apply processing, then use smooth windows on each segment before recombining. This ensures a smooth frequency transition between segments.
Homomorphic Signal Processing
•	The second nonlinear technique, homomorphic signal processing, deals with signals mixed through multiplication or convolution, not just addition. It transforms these nonlinear combinations into a linear form for easier separation.
•	Homomorphic signal processing separates nonlinearly combined signals by converting the problem to a linear one.
•	Nonlinear signal combinations, like multiplication and convolution, cannot be separated by linear filtering.
Homomorphic Processing of Multiplied Signals
•	For example, an AM radio signal's loudness may vary due to changing atmospheric conditions. This can be modeled as the audio signal a[n] multiplied by a slowly changing gain signal g[n]. While typically corrected by an automatic gain control (AGC) circuit, it can also be addressed using nonlinear DSP techniques.
•	Audio signals transmitted via AM radio waves, affected by changing atmospheric conditions, can be corrected with nonlinear DSP.
 
In homomorphic processing, signals combined by multiplication are turned into added signals using a logarithm. A conventional linear filter—typically a high-pass filter—removes the low-frequency gain component log(g[n]), leaving log(a[n]). Finally, the anti-logarithm (exponential) is applied to recover the desired signal a[n].

Homomorphic Processing of Convolved Signals
•	Homomorphic processing can separate signals mixed by convolution, such as removing echoes in audio. It uses a Fourier transform to turn convolution into multiplication, then a logarithm to convert multiplication into addition. After linear filtering, the inverse transform (anti-log + inverse FFT) restores the signal.
•	Interestingly, linear filtering here works in the frequency domain, similar to how signals are usually handled in the time domain—a reversal that introduces unique terms. For example, cepstrum (from “spectrum”) is the Fourier transform of the log of the spectrum. Filters become long-pass and short-pass, and terms like quefrency and liftering are used in this context.


 
 
•	The DFT contains a large spike around the “period” of the signal and some “low-frequency” components due to the amplitude modulation. 
•	A simple filter would then allow us to separate both components.
Homomorphic Signal Processing - problems
•	Homomorphic processing requires handling both negative and positive values using the complex logarithm.
•	Aliasing can occur when taking the logarithm, requiring significantly higher sampling rates.
•	Linearized signals may not be separable by linear filters if their spectra overlap, even if the original signals did not.
Homomorphic Signal Processing
•	Signals should be processed in a manner consistent with how they are formed.
•	Understanding how information is represented in the signals being processed is the first step in any DSP task.
