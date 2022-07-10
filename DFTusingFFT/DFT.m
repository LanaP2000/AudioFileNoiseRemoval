close all;
clear all;
clc;

n = [0:43];
x = 2 * cos((pi/2)*n) - 3 * sin((pi/11)*n);
dft_x = fft(x);

figure;
subplot(211);
plot(n, abs(dft_x))

X = myFFT(x);
subplot(212);
plot(n, abs(X))

function [my_fft] = myFFT(input_signal)
    N = size(input_signal, 2);
    even_signal = input_signal(1:2:end);
    odd_signal = input_signal(2:2:end);

    fft_even = fft(even_signal);
    fft_odd = fft(odd_signal);

    WN = exp(((-1i*2*pi)/N) * (0:N/2-1));

    my_fft = [fft_even + WN.*fft_odd fft_even - WN.*fft_odd];
end