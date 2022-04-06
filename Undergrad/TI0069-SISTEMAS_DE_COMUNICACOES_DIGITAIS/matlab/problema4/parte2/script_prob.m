%% Trabalho de SCD 
% ----------------------------
% Problema 4 - Probabilidade de Erro da Modulacao M-PSK
% script_prob.m
% 2021/03/27 - Lucas Abdalah
%

close all; clearvars; clc; % Clear the matlab ambient

%% Add necessary folder paths
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema2\'

%% General Setup
plot_figures = true;
save_figures = true;

%% Calculo da probabilidade do M-QAM
M = [4,8];                  % Simbolos das constelacoes
Es = [0.5,0.5];                  % Energia das constelacoes
Es_N0 = [0:2:20];               % Vetor de SNR em dB
N0 = Es'.*(10.^(-Es_N0/10));    % Ruido linear para cada constelacao

%% Calculo do erro teorico x Es/N0 em dB
erro_MPSK = zeros(size(M,2),length(Es_N0)); % Matriz onde cada linha recebe o erro em em funcao de Es/N0 em dB
for ii=1:size(M,2)
    for jj=1:length(Es_N0)
        erro_MPSK(ii,jj)= Pe_MPSK(M(ii),Es(ii),N0(ii,jj));
    end
end

%% Plotando erro teorico x Es/N0 em dB
if plot_figures == true
    h = figure;
    str = {'4PSK','8PSK'};
    semilogy_erro(Es_N0,erro_MPSK, str);
    title('$P(e)$ Teorico: $M$-PSK','interpreter','latex');
    ylim([1e-6,1e0]);
    grid on
    if save_figures == true
        saveas(h,['fig/Erro_Teorico_MPSK'],'pdf');
    end
end

% Salvando dados de SER e BER
filename = ['Problema4_M_PSK'];
save([filename,'.mat']);