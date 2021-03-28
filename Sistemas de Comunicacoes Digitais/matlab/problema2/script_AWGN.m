%% Trabalho de SCD 
% ----------------------------
% Problema 2 - Modulacao M-QAM em canal AWGN
% script_AWGN.m
% 2021/03/26 - Lucas Abdalah
%

% close all; clearvars; clc; % Clear the matlab ambient

% addpath 'problema1' % Local path
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab'
%% General Setup
close all; 
clearvars; 
clc

%% Calculo da probabilidade do M-QAM
M = [4,16,64];
Es = [1,5,21];
N_size=size(M,2);
Es_N0 = [0:2:20];  % em dB
N_dB = length(Es_N0);
N0 = Es'.*(10.^(-Es_N0/10));

erro_MQAM = zeros(N_size,N_dB);

for ii=1:N_size
    for jj=1:N_dB
        erro_MQAM(ii,jj)= Pe_MQAM(M(ii),Es(ii),N0(ii,jj));
    end
end

figure
str = {'4QAM','16QAM','64QAM'};
semilogy_erro(Es_N0,erro_MQAM, str);
title('$P(e)$ Teorico: $M$-QAM','interpreter','latex');
ylim([1e-6,1e0]);
grid on

% generate the noise term
% noise= sqrt(No/2)*(randn(N,1)+1i*randn(N,1));
% received signal with noise
% rx_signal= symb + noise; 

