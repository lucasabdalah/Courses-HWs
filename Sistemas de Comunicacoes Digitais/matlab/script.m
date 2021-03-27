%% Trabalho de SCD 
% Testes
% 2021/03/24 - Lucas Abdalah
%

close all; clearvars; clc; % Clear the matlab ambient

% addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema1'
addpath 'problema1' % Local path

M = 16;         % - Numero de simbolos da constelacao
L = 4096;       % - Tamanho da sequencia (bits)
K = log2(M);    % - Numero de bits/simbolo
N = L/K;        % - Numero de simbolos a serem transmitidos

% s_m = randi([0 1],L,1)';   % - Mensagem a transmitir

% Generating a gray const.
% G = gray_const(M,true);

x = 1 + 1j*1;
y = -1 + 0j*1;

d_E = sqrt( (real(x)-real(y))^2 + (imag(x)-imag(y))^2)



d = distancia_QAM(1);
M = 4
a = randn(1);
r = a + 1j*a
[symb,bits]=demapping_MQAM(r,M(1),d)
disp('---------------------------------------')
symb = mapping_MQAM(bits)
figure(hfig(1));
scatter(real(r),imag(r));
hold on 
hold off
