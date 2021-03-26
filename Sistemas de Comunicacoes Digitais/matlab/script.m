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

D = de2bi([0:M-1],log2(M),'left-msb');

G = zeros(M,log2(M));

% for ii=1:M
%     for jj=1:log2(M)        
%         G(ii,log2(M)) = D(ii,log2(M));
%     end
% end

% --------------------------------
bit_len = 4;
b = randi([0 1],bit_len,1)';
g = mybin2gray(b);
disp(['b ==> ', num2str(b)])
disp(['g ==> ', num2str(g)])

%dec2bin(str2num(sprintf('%d',x).'))
%
%
% Em_No = 30;     % - Razao sinal-ruido
%