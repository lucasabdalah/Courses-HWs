%% Trabalho de SCD 
% ----------------------------
% Problema 3 - Canal AWGN da Modulacao M-QAM
% script_prob.m
% 2021/03/27 - Lucas Abdalah
%

clearvars; close all; clc; % Clear the matlab ambient

%% Add necessary folder paths
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab'
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema1'

%% General Setup
plot_figures = true;
save_figures = true;

%% Calculo da probabilidade do M-QAM
M = 16;          % - Numero de simbolos da constelacao
K = log2(M);    % - Numero de bits/simbolo
L = 50*M*K;    % - Tamanho da sequencia (bits)
% L = 264000;     % - Tamanho da sequencia (bits)
N = L/K;        % - Numero de simbolos a serem transmiti

%% Calculo da energia media (E_m) e da distancia entre os simbolos
E_g = 1;                    % - Energia do pulso de transmissao
E_m = energia_MQAM(M,E_g);  % - Energia da constelacao MQAM
d = d_MQAM(M,E_m);          % - Distancia entre os simbolos
% coef = sqrt(2)*d; % Correction constant
%% Cadeia de bits
s_m = randi([0 1],L,1)';    % - Mensagem a transmitir

%% Simbolos a serem mapeados e demapeados
symb_map = zeros(1,N);      % - Vetor vazio para receber os simbolos mapeados

n=0;    % - Index for symbols iteration
for ii = 1:K:L
    n=n+1;   
    bits = s_m(ii:ii+K-1);
    symb_map(n) = mapping_MQAM(bits);
    % [symb_demap(n),y_m(ii:ii+K-1)]=demapping_MQAM(symb_map(n),M,d);
end 

%% Generate the noise term
% Es_N0=[0:2:20]; % em dB
Es_N0=25; % em dB
N0= E_m*10.^(-Es_N0/10);
noise= sqrt(N0/2)'*(randn(N,1)+1i*randn(N,1))';

%% Ma received signal with noise
rx_signal = symb_map + noise; 

if plot_figures == true
    h = figure;
    scatter_MQAM(const_MQAM(M,d),true,false,true);
    hold on;
    scatter_MQAM(rx_signal,false,true,false); 
    ax = max(abs(real(rx_signal)));
    ay = max(abs((imag(rx_signal))));
    xlim([-ax-d,ax+d]);
    ylim([-ay-d,ay+K*d/2]);
    title(['Constelacao ',num2str(M),'-QAM: ', num2str(Es_N0), 'dB de SNR']);
    str = {[num2str(M),'QAM constelation'],'Sinal Emitido'};
    legend(str,'Location','north', 'Orientation','Horizontal');
    if save_figures == true;
        saveas(h,['fig/',num2str(M),'QAM_',num2str(Es_N0),'dB'],'pdf');
    end
end

% pause

% %% Demapear os sinais
% symb_demap = zeros(1,N);    % - Vetor vazio para receber os simbolos demapeados
% y_m = zeros(1,L);           % - Vetor vazio para receber a mensagem

% n=0;    % - Index for symbols iteration
% for ii = 1:K:L
%     n=n+1;
%     [symb_demap(n),y_m(ii:ii+K-1)]=demapping_MQAM(rx_signal(n),M,d);
% end 

% %% Message length and Bit Error
% disp(['message length: ', num2str(L)]);
% disp(['bit error: ', num2str(sum(xor(s_m,y_m)))]);