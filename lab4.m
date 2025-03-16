clc; clear;
x = 0:127;
x = x / 1024;
y = -16*x + 3;
t = 0:127;
F = [y, zeros(1, 129), ones(1, 128), sin(2*pi*0.017*t) + 1];
N = 1024;
h = ifft(F, N, 'symmetric');
h_shifted = fftshift(h);

M_values = [40, 100, 300, 1000];

for i = 1:length(M_values)
    M = M_values(i);

    h_truncated = h_shifted(N/2 - M/2 : N/2 + M/2);
    window = hamming(M+1)';
    h_windowed = h_truncated .* window;

    H = fft(h_windowed, N);
    H_magnitude = abs(H(1:N/2+1));

    figure;
    title(['Filter Analysis for M = ', num2str(M)]);
    subplot(2,1,1);
    stem(h_windowed);
    title('Impulse Response');
    xlabel('Sample Index');
    ylabel('Amplitude');

    subplot(2,1,2);
    plot(linspace(0, 0.5, N/2+1), H_magnitude);
    title('Magnitude Response');
    xlabel('Normalized Frequency');
    ylabel('Magnitude');
end