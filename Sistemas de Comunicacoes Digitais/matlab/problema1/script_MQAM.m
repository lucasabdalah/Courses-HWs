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
plot_figures = false;
save_figures = false;

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
% Since the distances are equal, let's define a standard for all cases;
d=distancia_QAM(1);

%% --------------------------------------------------------------------
%% Plotando Constelacoes
if plot_figures == true
    for ii = 1:size(M,2)
        fprintf('>>Pressione qualquer tecla para ver a constelacao %1d-QAM \n', M(ii));
        pause;
        const_QAM = const_MQAM(M(ii),d);
        figure;
        h = scatter_MQAM(const_QAM,true);    
        % Generating Gray alphabet and constellation's plot.
        for jj = 1:M(ii)
            gray_alfabeto = gray_const(M(ii),false);
            str = {strjoin(string(gray_alfabeto(jj,:)))};
            text(real(const_QAM(jj)),imag(const_QAM(jj))+(ii/6),...
            str,...
            'FontSize', 6,...
            'HorizontalAlignment', 'center');
        end
        % 
        if save_figures == true
            % Uncomment to save the figures
            saveas(h,['fig/',num2str(M(ii)),'_QAM_plot'],'pdf');
        end
    end 
end 

pause

%% Creating the constellation for each M-QAM
const_4QAM = const_MQAM(M(1),d);
const_16QAM = const_MQAM(M(2),d);
const_64QAM = const_MQAM(M(3),d);

% const_16QAM
bits = [0,0,1,1];
symb = mapping_MQAM(bits)

%% --------------------------------------------------------------------
pause;

% Em_No = 30;     % - Razao sinal-ruido