%% Trabalho de SCD 
% ----------------------------
% Problema 2 - Probabilidade de Erro da Modulacao M-QAM
% script_prob.m
% 2021/03/26 - Lucas Abdalah
%

close all; clearvars; clc; % Clear the matlab ambient

%% Add necessary folder paths
% addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab'

%% General Setup
plot_figures = true;
save_figures = false;

%% Calculo da probabilidade do M-QAM
M = [4,16,64];                  % Simbolos das constelacoes
Es = [1,5,21];                  % Energia das constelacoes
Es_N0 = [0:2:20];               % Vetor de SNR em dB
% Es_N0 = [0:2:30];               % Vetor de SNR em dB para Erro_Teorico_MQAM_extended.pdf
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
        % saveas(h,['fig/Erro_Teorico_MQAM_extended'],'pdf');
    end
end

%% Salvando dados de SER e BER
filename = ['Problema2_M_QAM'];
save([filename,'.mat']);