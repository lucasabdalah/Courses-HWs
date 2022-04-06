%% General Parameters 
Fs = 1e3 ; % Hz
numSamples = 2*Fs; % Samples 
t = linspace(0,1,numSamples); % Time axis
%% Mixing Matrix
M = [0.1 0.2; 0.3 0.4]; % Mixing Matrix
%% Sources
w(1) = 2*pi*2; % Frequncy of the first signal
w(2) = 2*pi*20; % Frequency of the second signal
% S = ones(2,numSamples); % A matrix to receive the signals
S(1,:) = sin(w(1)*t);
S(2,:) = sin(w(2)*t);
% Plot the sources 
figure;
subplot(2,1,1);
plot(t,S(1,:)); % plot
title('Sources')
ylabel('Amplitude')
grid on;
subplot(2,1,2); 
plot(t,S(2,:));
xlabel('Time (s)')
grid on;
%% Mixed Signals
X = M * S; % Mixing using M
% Plot the signals
figure
subplot(2,1,1);
plot(t,X(1,:)); % plot
title('Mixed Signals')
ylabel('Amplitude')
grid on;
subplot(2,1,2); 
plot(t,X(2,:));
xlabel('Time (s)')
grid on;
%% Unmixing Signals
Y = M \ X; % Unmixing using M
% Plot the signals
figure
subplot(2,1,1);
plot(t,Y(1,:)); % plot
title('Source Estimation')
ylabel('Amplitude')
grid on;
subplot(2,1,2); 
plot(t,Y(2,:));
xlabel('Time (s)')
grid on;
