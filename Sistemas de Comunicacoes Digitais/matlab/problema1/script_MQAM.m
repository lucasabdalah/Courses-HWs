%% Trabalho de SCD 
% ----------------------------
% Problema 1 - Modulacao M-QAM 
% script_MQAM.m
% 2021/03/24 - Lucas Abdalah
%

close all; clearvars; clc; % Clear the matlab ambient

% addpath 'problema1' % Local path
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab'
%% General Setup
% save_figures = true;
%% Color definitions for graphics
% yellow = [0.9290 0.6940 0.1250];
% black  = [0 0 0];

%% Problema 1 - Calculo da energia media (E_m) e da distancia entre os simbolos
M = [4, 16, 64];    % - Numero de simbolos da constelacao
E_g = 1;            % - Energia do pulso de transmissao

energia_media_QAM = zeros(1,size(M,2));
distancia_QAM = zeros(1,size(M,2));

for ii = 1:size(M,2)
    energia_media_QAM(1,ii) = energia_MQAM(M(ii),E_g);
    distancia_QAM(1,ii) = d_MQAM(M(ii),energia_media_QAM(1,ii));
    fprintf('--------------------- Modulacao %1d-QAM --------------------- \n', M(ii));
    fprintf('Energia (E_media) = %1d \n', energia_media_QAM(1,ii));    
    fprintf('Distancia (d) = %1d \n\n', distancia_QAM(1,ii));
end
fprintf('------------------------------------------------------------\n');

%% --------------------------------------------------------------------
%% Creating the constellation for each M-QAM
const_4QAM = const_MQAM(M(1),distancia_QAM(1));
const_16QAM = const_MQAM(M(2),distancia_QAM(2));
const_64QAM = const_MQAM(M(3),distancia_QAM(3));

figure
pre_gray_16QAM = const_16QAM;
h = scatter_QAM(pre_gray_16QAM,true);
ii=2;
for jj = 1:M(ii)
    % bin_alfabeto = de2bi(0:M(ii)-1,'left-msb');
    bin_alfabeto = gray_const(M(2),true);
    str = {strjoin(string(bin_alfabeto(jj,:)),'\n')}; % Create
    text(real(pre_gray_16QAM(jj)),imag(pre_gray_16QAM(jj))+0.3,...
    str);
end

pause 
%% Plotando Constelacoes
if save_figures == true
    
    for ii = 1:size(M,2)
        fprintf('>>Pressione qualquer tecla para ver a constelacao %1d-QAM \n', M(ii));
        pause;
        const_QAM = const_MQAM(M(ii),distancia_QAM(ii));
        figure;
        h = scatter_QAM(const_QAM,true);
        % Uncomment to save the figures
        saveas(h,['fig/',num2str(M(ii)),'_QAM_plot'],'pdf');
        %
        % This part is to display the symbol's bits
        % To make this part work is necessary to organize the M-QAM
        % modulation to sort the binary constellation and convert it into encoding
        % for jj = 1:M(ii)
        %     bin_alfabeto = de2bi(0:M(ii)-1,'left-msb');
        %     str = {strjoin(string(bin_alfabeto(jj,:)),'\n')}; % Create
        %     text(real(const_QAM(jj)),imag(const_QAM(jj))+0.3,...
        %     str);
        % end
    end 
end 


%% --------------------------------------------------------------------
pause;
% close all;

%
% scatter(real(vec(const_4QAM)), imag(vec(const_4QAM)))
% A construcao da constelacao segue a logica da esquerda para direita (->) debaixo pra cima (^)
% Em_No = 30;     % - Razao sinal-ruido