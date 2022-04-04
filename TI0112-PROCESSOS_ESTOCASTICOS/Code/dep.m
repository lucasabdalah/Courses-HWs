% Sine signal
Fs = 1000;
t = 0:1/Fs:1-1/Fs;
x_sine = sin(2*pi*150*t) + sin(2*pi*180*t);

N = length(x_sine);
xdft = fft(x_sine);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x_sine):Fs/2;

% Random signal
%X = exprnd(2,1,N);
X = randn(1,N);
r_X = xcorr(X)/N;
S_r = fft(r_X);
S_r = S_r(1:N/2+1);
S_r = (1/(Fs*N)) * abs(S_r).^2;
S_r(2:end-1) = 2*S_r(2:end-1);


% Sine + noise signal
x_noisy = x_sine + 0.5*X;
S_x = xcorr(x_noisy)/N;
x_noisy_dft = fft(x_noisy);
x_noisy_dft = x_noisy_dft(1:N/2+1);
psdx_noisy = (1/(Fs*N)) * abs(x_noisy_dft).^2;
psdx_noisy(2:end-1) = 2*psdx_noisy(2:end-1);
freq = 0:Fs/length(x_noisy):Fs/2;




figure
plot(freq,10*log10(psdx));
hold on
plot(freq,10*log10(S_r),'r');
plot(freq,10*log10(psdx_noisy),'k');
legend('Seno','Aleatorio','Sinal+ruido');

