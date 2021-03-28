%% Trabalho de SCD 
% ----------------------------
% Problema 2 - Modulacao M-QAM em canal AWGN
% script_AWGN.m
% 2021/03/26 - Lucas Abdalah
%

close all; clearvars; clc; % Clear the matlab ambient

addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab'
% addpath ''C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema1'
%% General Setup
plot_figures = false;
save_figures = false;

%% Calculo da probabilidade do M-QAM
M = [4,16,64];                  % Simbolos das constelacoes
Es = [1,5,21];                  % Energia das constelacoes
Es_N0 = [0:2:20];               % Vetor de SNR em dB
N0 = Es'.*(10.^(-Es_N0/10));    % Ruido linear para cada constelacao

%% Calculo do erro teorico x Es/N0 em dB
erro_MQAM = zeros(size(M,2),length(Es_N0)); % Matriz onde cada linha recebe o erro em em funcao de Es/N0 em dB
for ii=1:size(M,2)
    for jj=1:length(Es_N0)
        erro_MQAM(ii,jj)= Pe_MQAM(M(ii),Es(ii),N0(ii,jj));
    end
end

%% Plotando erro teorico x Es/N0 em dB
if plot_figures == true
    h = figure;
    str = {'4QAM','16QAM','64QAM'};
    semilogy_erro(Es_N0,erro_MQAM, str);
    title('$P(e)$ Teorico: $M$-QAM','interpreter','latex');
    ylim([1e-6,1e0]);
    grid on
    if save_figures == true
        saveas(h,['fig/Erro_Teorico_MQAM'],'pdf');
    end
end


%
% generate the noise term
% noise= sqrt(No/2)*(randn(N,1)+1i*randn(N,1));
% received signal with noise
% rx_signal= symb + noise; 