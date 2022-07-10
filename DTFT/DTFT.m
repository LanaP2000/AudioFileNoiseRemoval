close all
clear all
clc

%% DTFT of x[n] = (0.5)^n(u[n])
omega1 = 0:pi/1024:pi-(pi/1024);
b = [1];
a = [1 -0.5];

[X, omega1] = freqz(b, a, omega1);

figure; subplot(2, 1, 1)
plot(omega1, abs(X))
subplot(2, 1, 2)
plot(omega1, unwrap(angle(X)))

%% y[n] + 0.6y[n -1] + 0.03y[n - 2] - 0.01y[n - 3] = x[n] + 2x[n -1] - x[n - 2]
n = 0:10;
% Impulse Response
b = [1 2 -1];
a = [1 0.6 0.03 -0.01];

hnz = filter(b, a, [1, zeros(1, length(n) - 1)]);

figure;
stem(n, hnz, 'filled')
xlabel('Time index n')
ylabel('h[n]')

%% 
omega2 = 0:pi/10:pi-(pi/10)
b = [1 -2 -1];
a = [1 0.1 -0.025 0.05];
ynz = filter(b, a, [1, ones(1, length(omega2)-1)]);

figure;
subplot(311)
stem(omega2, ynz, 'filled')

[aa bb] = invfreqz(ynz,omega2,2,4)
[h,omega2] = freqz(aa,bb,omega2);

subplot(312)
stem(omega2, h, 'filled')

err = h - ynz
subplot(313)
stem(omega2, err, 'filled')
