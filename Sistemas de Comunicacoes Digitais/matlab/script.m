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

%% Problema 1 - Calculo da energia media (E_media)
E_g = 1;    % - Energia do pulso de transmissao
for ii = 1:3
    fprintf('--------------------- Modulacao %1d-QAM --------------------- \n', 4^ii);
    fprintf('Energia (E_media) = %1d \n\n', energia_MQAM(4^ii,E_g));
end

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
% scatter()
%
% figure
% hold on
% for i = 1:4:16
%     display(['i=', num2str(i)])
%     slice = y(i:i+3);
%     display(num2str(slice));
%     scatter(real(slice), imag(slice),'filled');
%     pause;
% end
% legend
% xlim([-4,4]);
% ylim([-4,4]);
% hold off

%     text(real(slice(i)), imag(slice(i)),{'\leftarrow a'});
% scatterplot(y)
% https://www.mathworks.com/matlabcentral/fileexchange/14809-m-qam-modulation-and-demodulation

% scatter(y,'AlignVertexCenters','on','MarkerSize',25,'Marker','.',...
%     'Color',[0 0 1]);