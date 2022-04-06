%% AP2 de PES
% Questao 2
%
% Q2_script.m
%
% 2021/08/23 - Lucas Abdalah

close all; clearvars; clc; % Clear the matlab ambient

% To reproduce the same results
rng('default');
pausetime = 6;

%% Sequencia de 1000 amostras N(0,1)
N = 1e3;

%% Sequencia com 2N amostras 
% N = 2*N; %<-- Desfazer comentario para resultados do item a)
% N = 20*N; %<-- Desfazer comentario para resultados do item d)

H_0 = randn(1,N);


%% Propriedades Estatisticas
    fprintf('Propriedades Estatisticas de H_0 com N = %1d \n-------------------------------------------\n', N);
    fprintf('Media = %1d \n', mean(H_0));
    fprintf('Variancia = %1d \n', var(H_0));
    fprintf('-------------------------------------------\n');

pause(pausetime);
clc;
%% item a) Limiar de Bayes \eta_{B}
    % Teste de H_0
    fprintf('--> item a) Limiar de Bayes eta_{B} <--\n-------------------------------------------\n');
    limiter_B = -1.0404; 
    detector_B = 0;

    for ii = 1:1:N
        if H_0(ii) >= limiter_B
            detector_B = detector_B + 1;
        end
    end

    P_10_B = (detector_B/N);
    fprintf('Teste de H_0 com %1d amostras:\n', N);
    fprintf('P_10 = %2.4f \n', P_10_B);
    fprintf('-------------------------------------------\n');
    % Calculo numerico da integral
    P_10_B_int = qfunc(limiter_B);
    fprintf('Calculo numerico da integral (Limiar->inf):\n');
    fprintf('P_10 = %2.4f \n-------------------------------------------\n', P_10_B_int);
    fprintf('Erro Percentual = %2.4f \n-------------------------------------------\n', 100*abs(((P_10_B - P_10_B_int)/P_10_B_int)));

pause(pausetime);
clc;
%% item b) Limiar MAP \eta_{MAP}
    % Teste de H_0
    fprintf('--> item b) Limiar MAP eta_{MAP} <--\n-------------------------------------------\n');
    limiter_MAP = -0.3473; 
    detector_MAP = 0;
    
    for ii = 1:1:N
        if H_0(ii) >= limiter_MAP
            detector_MAP = detector_MAP + 1;
        end
    end
    
    P_10_MAP = (detector_MAP/N);
    fprintf('Teste de H_0 com %1d amostras:\n', N);
    fprintf('P_10 = %2.4f \n', P_10_MAP);
    fprintf('-------------------------------------------\n');
    % Calculo numerico da integral
    P_10_MAP_int = qfunc(limiter_MAP);
    fprintf('Calculo numerico da integral (Limiar->inf):\n');
    fprintf('P_10 = %2.4f \n-------------------------------------------\n', P_10_MAP_int);
    fprintf('Erro Percentual = %2.4f \n-------------------------------------------\n', 100*abs(((P_10_MAP - P_10_MAP_int)/P_10_MAP_int)));

pause(pausetime);
clc;
%% item c) Limiar Minimax \eta_{MN}
    % Teste de H_0
    fprintf('--> item c) Limiar Minimax eta_{MN} <--\n-------------------------------------------\n');
    limiter_MN = 0.5; 
    detector_MN = 0;
    
    for ii = 1:1:N
        if H_0(ii) >= limiter_MN
            detector_MN = detector_MN + 1;
        end
    end

    P_10_MN = (detector_MN/N);
    fprintf('Teste de H_0 com %1d amostras:\n', N);
    fprintf('P_10 = %2.4f \n', P_10_MN);
    fprintf('-------------------------------------------\n');
    % Calculo numerico da integral
    P_10_MN_int = qfunc(limiter_MN);
    fprintf('Calculo numerico da integral (Limiar->inf):\n');
    fprintf('P_10 = %2.4f \n-------------------------------------------\n', P_10_MN_int);
    fprintf('Erro Percentual = %2.4f \n-------------------------------------------\n', 100*abs(((P_10_MN - P_10_MN_int)/P_10_MN_int)));
pause(pausetime);
clc;
%% item d) Limiar Neyman-Pearson \eta_{NP}
    % Teste de H_0
    fprintf('--> item d) Limiar Neyman-Pearson eta_{NP} <--\n-------------------------------------------\n');
    limiter_NP = 2.5; 
    detector_NP = 0;
    
    for ii = 1:1:N
        if H_0(ii) >= limiter_NP
            detector_NP = detector_NP + 1;
        end
    end

    P_10_NP = (detector_NP/N);
    fprintf('Teste de H_0 com %1d amostras:\n', N);
    fprintf('P_10 = %2.4f \n', P_10_NP);
    fprintf('-------------------------------------------\n');
    % Calculo numerico da integral
    P_10_NP_int = qfunc(limiter_NP);
    fprintf('Calculo numerico da integral (Limiar->inf):\n');
    fprintf('P_10 = %2.4f \n-------------------------------------------\n', P_10_NP_int);
    fprintf('Erro Percentual = %2.4f \n-------------------------------------------\n', 100*abs(((P_10_NP - P_10_NP_int)/P_10_NP_int)));
        
    