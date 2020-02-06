function [y,fs] =  trompete(filename)
    clc; % Limpa a janela de comandos
    
    % Recebendo o audio com 2 canais
    [file, fs] = audioread(filename); % Le o audio contido no mesmo diretorio - y carrega os dados e fs a freq de amostragem
    t = linspace(0,length(file)/fs,length(file)); % Tempo do Audio
    
    % Grafico 1 - Amostras x Tempo audio stereo
    figure(1);
	('GRAFICO 1 - AMOSTRAS x TEMPO');
    % Plot do Sinal contra o Tempo
	plot(t, file);
    title ("Sinal do Audio (fs = 44.1KHz)"); %title of this graphic
	xlabel ("tempo (s)"); 
    ylabel ("Amplitude (Normalizada)");
    %x is the time index in seconds, y is the amplitude's samples 
	grid on; %turn on the grid
    
    % ESQUERDA
    E = file(:,1); % Canal Esquerdo
    n = length(E);
    % FFT para analise
    [p, q] = fftedit(E, fs);
    
    % Subplot do AUDIO NO TEMPO E NA FREQUENCIA
    figure(2)
    subplot(2,1,1); 
    plot(t, E);  % Esquerda contra o Tempo
    title ("Subplot 1: Esquerda (fs = 44.1KHz)"); %title of this graphic
	xlabel ("tempo (s)"); % x is the time index in seconds,
    ylabel ("Amplitude (Normalizada)"); % y is the amplitude's samples 
    grid on; %turn on the grid
    subplot(2,1,2);
    plot(q,p); % FFT em Hertz
    title("Subplot 2: Frequencia em Hertz");
    ylabel("Amplitude (Não Normalizada)")
    grid on; %turn on the grid
    
    % Esquerda
    figure(3);
    spectrogram(E, fs); % Analise na Frequencia
       
    % 0RDEM 10
    delta1 = 1e-3; % Limite da Banda de Passagem 
    deltas = 1e-3; % Limite da Banda de Rejeicao
    wp = 0.218*pi; % Corte 1
    ws = 0.454*pi; % Corte 2 
    
    % Calculo da ordem N e do omegac
    auxdelta = log((1/deltas)^2 - (1/1-delta1)^2); % mudanca de base
    ausw = log(ws/wp); % % mudanca de base
    disp('Ordem N:');
    ORDEM10 = 0.5*(auxdelta/ausw); %
    disp('Ordem N inteiro:'); disp(ORDEM10);% 
    ORDEM10 = ceil(ORDEM10); disp(ORDEM10);% 

    omegac = (wp/(((1/(1-delta1))^2 - 1)^0.5*ORDEM10));
    disp('Omegac:'); disp(omegac);
    
	% FILTRO
    fc = 4800; % Frequencia de corte calculada
    [a,b] = butter(ORDEM10,fc/(fs/2)); % Filtro ButterWorth
    yE = filter(a,b,E); %audio filtrado
    % sound(yE,fs);
    
    % PLOT DO FILTRO
    figure(4);
    freqz(a,b,fs);
    title('Filtro Passa Baixas Projetado');
    
    % POLOS NO PLANO Z
    figure (5);
    Hz = filt(a,b,1);
    pzmap(Hz);
    
    % POLOS NO PLANO S
    figure (6)
    Hs = d2c(Hz);
    pzmap(Hs);
    
    % Direita
    D = file(:,2); % Canal Direito
    
    % FFT para analise
    [p, q] = fftedit(D, fs);
    
    % Subplot do AUDIO NO TEMPO E NA FREQUENCIA
    figure(7)
    subplot(2,1,1); 
    plot(t, D);  % Direita contra o Tempo
    title ("Subplot 1: Direita (fs = 44.1KHz)"); %title of this graphic
	xlabel ("tempo (s)"); % x is the time index in seconds,
    ylabel ("Amplitude (Normalizada)"); % y is the amplitude's samples 
    grid on; %turn on the grid
    subplot(2,1,2);
    plot(q,p); % FFT em Hertz
    title("Subplot 2: Frequencia em Hertz");
    ylabel("Amplitude (Não Normalizada)")
    grid on; %turn on the grid
    
    figure(8);
    spectrogram(D, fs); % Analise na Frequencia

    yD = filter(a,b,D); %audio filtrado
    
    y = [yE yD];
    
    sound(y,fs);
    audiowrite('q1fitro.wav',y,fs);
    
	% 0RDEM 5
    delta1 = 1e-1; % Limite da Banda de Passagem 
    deltas = 1e-1; % Limite da Banda de Rejeicao
    wp = 0.218*pi; % Corte 1
    ws = 0.480*pi; % Corte 2 
    % Calculo da ordem N e do omegac
    auxdelta = log((1/deltas)^2 - (1/1-delta1)^2); % mudanca de base
    ausw = log(ws/wp); % % mudanca de base
    disp('Ordem N:');
    ORDEM5 = 0.5*(auxdelta/ausw); %
    disp('Ordem N inteiro:'); disp(ORDEM5);% 
    ORDEM5 = ceil(ORDEM5); disp(ORDEM5); 
 
    % FILTRO
    fc = 4800; % Frequencia de corte calculada
    [c,d] = butter(ORDEM5,fc/(fs/2)); % Filtro ButterWorth
    xE = filter(c,d,E); %audio filtrado
    xD = filter(c,d,D); %audio filtrado
    
    figure(9);
    freqz(a,b,fs);
    hold on;
    freqz(c,d,fs);
    title('Filtro Passa Baixas Projetado');
    
    x = [xE xD];
     
    audiowrite('q1fitroruim.wav',x,fs);
    
    
    
    
    
    [p, q] = fftedit(yE, fs);
    % Subplot do AUDIO NO TEMPO E NA FREQUENCIA
    figure(10)
    subplot(2,1,1); 
    plot(t, yE);  % Esquerda contra o Tempo
    title ("Subplot 1: Esquerda (fs = 44.1KHz)"); %title of this graphic
	xlabel ("tempo (s)"); % x is the time index in seconds,
    ylabel ("Amplitude (Normalizada)"); % y is the amplitude's samples 
    grid on; %turn on the grid
    subplot(2,1,2);
    plot(q,p); % FFT em Hertz
    title("Subplot 2: Frequencia em Hertz");
    ylabel("Amplitude (Não Normalizada)")
    grid on; %turn on the grid
    
    
    [p, q] = fftedit(yD, fs);
    % Subplot do AUDIO NO TEMPO E NA FREQUENCIA
    figure(11)
    subplot(2,1,1); 
    plot(t, yD);  % Direita contra o Tempo
    title ("Subplot 1: Direita (fs = 44.1KHz)"); %title of this graphic
	xlabel ("tempo (s)"); % x is the time index in seconds,
    ylabel ("Amplitude (Normalizada)"); % y is the amplitude's samples 
    grid on; %turn on the grid
    subplot(2,1,2);
    plot(q,p); % FFT em Hertz
    title("Subplot 2: Frequencia em Hertz");
    ylabel("Amplitude (Não Normalizada)")
    grid on; %turn on the grid
    
    
    
    
    
end