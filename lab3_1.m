%task1
clc; clear; close all;
N = 101;  
x = -50:50;
ft_values = [0.1, 0.25, 0.4];

figure;
hold on; grid on;
title('Impulse Response of Low-Pass Filter for Different Cutoff Frequencies');
xlabel('Samples'); ylabel('Amplitude');

for i = 1:length(ft_values)
    ft = ft_values(i);
    h = sinc(2 * ft * x);

    window = hamming(N)';
    h_windowed = h .* window;
    h_windowed = h_windowed / sum(h_windowed);
    
    plot(x, h_windowed, 'DisplayName', sprintf('ft = %.2f', ft));
end

legend show;
hold off;