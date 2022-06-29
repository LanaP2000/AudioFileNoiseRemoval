
[audio_raw,Fs]=audioread('Star_Wars.wav'); % Audio File Input
[N, P] = size(audio_raw); % File Size checking

% Plot in Time Domain
ts = 1/Fs; % Time Calculations for plotting
tmax=(N-1)*ts;
t=0:ts:tmax;

figure(1); % Figure 1 to compare input and output
subplot(211);
plot(t,audio_raw) % Plotting Raw Signal
title('Orignal Signal (Time Domain)'); % Labeling
xlabel('Time (t)');
ylabel('Amplitude');

% Fourier Transform
f=-Fs/2:Fs/(N-1):Fs/2;
z=fftshift(fft(audio_raw)); % Plotting FFT Spectrum of Signal

figure(2);
subplot(221);
plot(f,abs(z))
title('Fourier Transform');
xlabel('Frequency');
ylabel('Magnitude');

%% First Filter Design
N = 4000; % Order
Fc1 = 5000; % First Cutoff Frequency
Fc2 = 5300; % Second Cutoff Frequency
flag = 'scale'; % Sampling Flag
Beta = 0.5; % Window Parameter

% Create the window vector for the design algorithm
win = kaiser(N+1, Beta);

% Calculate the coefficients using the FIR1 function
b = fir1(N, [Fc1 Fc2]/(Fs/2), 'stop', win, flag);
Hd = dfilt.dffir(b);
freqz(Hd,N) % Frequency Response of First Filter

% Applying Filter to the signal
DataOut = filter(Hd,audio_raw);
z=fftshift(fft(DataOut));

figure(2);
subplot(222);
plot(f,abs(z)) % Plotting first output
title('Spike Removal at 5134 Hz with BandStop Filter');
xlabel('Frequency');
ylabel('Magnitude');

%% Second Filter Design
N = 1000; % Order
Fc1 = 3000; % First Cutoff Frequency
Fc2 = 4000; % Second Cutoff Frequency
flag = 'scale'; % Sampling Flag
Beta = 0.5; % Window Parameter

% Create the window vector for the design algorithm
win = kaiser(N+1, Beta);

% Calculate the coefficients using the FIR1 function
b = fir1(N, [Fc1 Fc2]/(Fs/2), 'stop', win, flag);
Hd = dfilt.dffir(b);
freqz(Hd,N) % Frequency Response of Second Filter
DataOut1 = filter(Hd,DataOut); % Applying Second Filter to Signal
z=fftshift(fft(DataOut1));

figure(2);
subplot(223) % Plotting Results
plot(f,abs(z))
title('Noise Removal @ 3000-4000 Hz');
xlabel('Frequency');
ylabel('Magnitude');

%% Third Filter Design
N = 5000; % Order
Fc = 200; % Cutoff Frequency
flag = 'scale'; % Sampling Flag
Beta = 0.5; % Window Parameter

% Create the window vector for the design algorithm
win = kaiser(N+1, Beta);
% Calculate the coefficients using the FIR1 function
b = fir1(N, Fc/(Fs/2), 'high', win, flag);
Hd = dfilt.dffir(b);
freqz(Hd,N) % Frequency Response of Third Filter
DataOut2 = filter(Hd,DataOut1); % Applying filter eliminate DC Component
z=fftshift(fft(DataOut2));

figure(2);
subplot(224); % Plotting Final Output
plot(f,abs(z))
title('Low Frequency Filter');
xlabel('Frequency');
ylabel('Magnitude');

figure(1);
subplot(212);
plot(t,DataOut2) % Plotting Output in Time Domain
title('Filtered Signal (Time Domain)');
xlabel('Time (t)');
ylabel('Amplitude');

audiowrite('Clean_File.wav',DataOut2,Fs)
P = audioplayer(DataOut2,Fs); % Playback the audio
play(P)


