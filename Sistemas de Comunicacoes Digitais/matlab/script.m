%% Trabalho de SCD 
% Problema 1 - Modulacao M-QAM 
% 2021/03/24 - Lucas Abdalah
%

close all; clearvars; clc; % Clear the matlab ambient

% addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema1'
addpath 'problema1' % Local path

M = 16;         % - Numero de simbolos da constelacao
L = 4096;       % - Tamanho da sequencia (bits)
K = log2(M);    % - Numero de bits/simbolo
N = L/K;        % - Numero de simbolos a serem transmitidos

s_m = randi([0 1],L,1)';   % - Mensagem a transmitir

% Generating a gray const.
G = gray_const(M,true);

% vec2flip = 1:M;
% Hi, I have a matrix xx which is a 6x20 matrix, and I want to flip only the even rows. Im trying this but is doesn't work:

% xx = reshape([1:12],[3,4])
% xx_flip = xx;
% xx_flip(2:2:end,:) = fliplr(xx(2:2:end,:))

% Funcao para flipar linha pares da matriz e mapear no codigo de Gray
% prepare_const2gray(const)

sqrt(M) % Cada linha apresenta sqrt(M) simbolos, entao para ordenar de modo a ter a codificao de Gray, teremos que inverter as linhas pares da constelacao

