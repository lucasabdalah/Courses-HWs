%% Trabalho de SCD 
% Problema 1 - Modulacao M-QAM 
% 2021/03/24 - Lucas Abdalah
%

close all; clearvars; clc; % Clear the matlab ambient

% addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema1'
addpath 'problema1' % Local path

M = 64;         % - Numero de simbolos da constelacao
L = 4096;       % - Tamanho da sequencia (bits)
K = log2(M);    % - Numero de bits/simbolo
N = L/K;        % - Numero de simbolos a serem transmitidos

s_m = randi([0 1],L,1)';   % - Mensagem a transmitir



%% Grey Code
% function b = gray2bin(g)
%     b(1) = g(1);
%     for i = 2 : length(g);
%         x = xor(str2num(b(i-1)), str2num(g(i)));
%         b(i) = num2str(x);
%     end
%
% Em_No = 30;     % - Razao sinal-ruido
%
%

% https://www.mathworks.com/matlabcentral/fileexchange/14809-m-qam-modulation-and-demodulation

% scatter(y,'AlignVertexCenters','on','MarkerSize',25,'Marker','.',...
%     'Color',[0 0 1]);