%% Trabalho de SCD 
% Testes
% 2021/03/24 - Lucas Abdalah
%

close all; clearvars; clc; % Clear the matlab ambient

% addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema1'
addpath 'problema1' % Local path
addpath 'problema4' % Local path
addpath 'problema4\parte1' % Local path

M = 4;         % - Numero de simbolos da constelacao
% L = 4096;       % - Tamanho da sequencia (bits)
K = log2(M);    % - Numero de bits/simbolo
% N = L/K;        % - Numero de simbolos a serem transmitidos

% s_m = randi([0 1],L,1)';   % - Mensagem a transmitir

% Generating a gray const.
% G = gray_const(M,true);

%% Problema 5 - Calculo da energia media (E_m) e da distancia entre os simbolos