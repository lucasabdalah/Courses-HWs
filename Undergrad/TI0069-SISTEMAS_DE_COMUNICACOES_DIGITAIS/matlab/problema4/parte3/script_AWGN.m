%% Trabalho de SCD 
% ----------------------------
% Problema 4 - Canal AWGN da Modulacao M-PSK
% script_AWGN.m
% 2021/03/27 - Lucas Abdalah
%

clearvars; close all; clc; % Clear the matlab ambient

%% Add necessary folder paths
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab'
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema1'
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema2'

%% General Setup
plot_figures = false;
save_figures = false;

%% Calculo da probabilidade do M-PSK
M = 8;         % - Numero de simbolos da constelacao
K = log2(M);    % - Numero de bits/simbolo
% L = 50*M*K;    % - Tamanho da sequencia (bits)
L = 264000;     % - Tamanho da sequencia (bits)
N = L/K;        % - Numero de simbolos a serem transmiti

%% Calculo da energia media (E_s) e da distancia entre os simbolos
E_g = 1;                    % - Energia do pulso de transmissao
E_s = energia_MPSK(M,E_g);  % - Energia da constelacao MPSK
d = d_MPSK(M,E_s);          % - Distancia entre os simbolos

%% Cadeia de bits
s_m = randi([0 1],L,1)';    % - Mensagem a transmitir

%% Simbolos a serem mapeados e demapeados
symb_map = zeros(1,N);      % - Vetor vazio para receber os simbolos mapeados

n=0;    % - Index for symbols iteration
for ii = 1:K:L
    n=n+1;   
    bits = s_m(ii:ii+K-1);
    symb_map(n) = mapping_MPSK(bits);
    % [symb_demap(n),y_m(ii:ii+K-1)]=demapping_MPSK(symb_map(n),M,d);
end 

if plot_figures == true

    %% Generate the noise term
    Es_N0=25; % em dB
    N0= E_s*10.^(-Es_N0/10);
    noise= sqrt(N0/2)'*(randn(N,1)+1i*randn(N,1))';

    %% Ma received signal with noise
    rx_signal = symb_map + noise;

    h = figure;
    scatter_MPSK(const_MPSK(M,E_s),true,false,true);
    hold on;
    scatter_MPSK(rx_signal,false,true,false); 
    ax = max(abs(real(rx_signal)));
    ay = max(abs((imag(rx_signal))));
    xlim(1.3*[-ax,ax]);
    ylim(1.3*[-ay,ay]);
    title(['Constelacao ',num2str(M),'-PSK: ', num2str(Es_N0), 'dB de SNR']);
    str = {[num2str(M),'PSK constelation'],'Sinal Emitido'};
    legend(str,'Location','north', 'Orientation','Horizontal');
    if save_figures == true;
        saveas(h,['fig/',num2str(M),'PSK_',num2str(Es_N0),'dB'],'pdf');
    end
    
    %% Demapear os sinais
    symb_demap = zeros(1,N);    % - Vetor vazio para receber os simbolos demapeados
    y_m = zeros(1,L);           % - Vetor vazio para receber a mensagem
    
    n=0;    % - Index for symbols iteration
    for ii = 1:K:L
        n=n+1;
        [symb_demap(n),y_m(ii:ii+K-1)]=demapping_MPSK(rx_signal(n),M,E_s);
    end     
end

pause

%% SER and BER plot -----------------
Es_N0=[0:2:20]; % em dB

% Vetores que recebem os valores de SER e BER
SER_curve = zeros(1,length(Es_N0));
BER_curve = zeros(1,length(Es_N0));

for ii = 1:length(Es_N0);
    %% Noise generation
    disp(['it:',num2str(ii)])
    N0= E_s*10^(-Es_N0(ii)/10);
    noise= sqrt(N0/2)'*(randn(N,1)+1i*randn(N,1))';

    %% Mapped signal with noise
    rx_signal = symb_map + noise;

    %% Demapear os sinais
    symb_demap = zeros(1,N);    % - Vetor vazio para receber os simbolos demapeados
    y_m = zeros(1,L);           % - Vetor vazio para receber a mensagem

    n=0;    % - Index for symbols iteration
    for jj = 1:K:L
        n=n+1;
        [symb_demap(n),y_m(jj:jj+K-1)]=demapping_MPSK(rx_signal(n),M,E_s);
    end

    %% parte
    BER = sum(xor(s_m,y_m))/L;
    diff_vec=symb_map-symb_demap;
    SER= sum(diff_vec~=0)/N;
    disp(['message length: ', num2str(L)]);
    disp(['bit error: ', num2str(sum(xor(s_m,y_m)))]);
    disp(['BER: ', num2str(BER)]);
    disp(['SER: ', num2str(SER)]);
    %% parte 
    BER_curve(1,ii) = BER;
    SER_curve(1,ii) = SER;
end

%% Gerando a curva SER
h = figure;
str = {[num2str(M),'PSK']};
semilogy_erro(Es_N0,SER_curve,str);
ylabel('SER');
title(['SER ', num2str(M),'-PSK'],'interpreter','latex');
saveas(h,['fig/',num2str(M),'PSK_','SER'],'pdf');

%% Gerando a curva BER
h = figure;
str = {[num2str(M),'PSK']};
semilogy_erro(Es_N0,BER_curve,str);
ylabel('BER');
title(['BER ', num2str(M),'-PSK'],'interpreter','latex');
saveas(h,['fig/',num2str(M),'PSK_','BER'],'pdf');

pause 
close all;

%% Salvando dados de SER e BER
filename = ['Problema4_',num2str(M),'_PSK'];
save([filename,'.mat']);