%% Trabalho de SCD 
% Problema 1 - Modulacao M-QAM 
% 2021/03/24 - Lucas Abdalah
%

close all; clearvars; clc; % Clear the matlab ambient

% addpath 'problema1' % Local path

%% Color definitions for graphics
yellow = [0.9290 0.6940 0.1250];
black  = [0 0 0];

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

%% Plotando Constelacoes
for ii = 1:size(M,2)
    fprintf('>>Pressione qualquer tecla para ver a constelacao %1d-QAM \n', M(ii));
    pause;
    const_QAM = const_MQAM(M(ii),distancia_QAM(ii));
    figure;
    h = scatter(real(const_QAM),imag(const_QAM), 'filled',...
    'MarkerEdgeColor',yellow,... 
    'MarkerFaceColor',yellow,...
    'LineWidth',1.5);
    % h.SizeData = 200;
    xlabel('In-Phase');
    ylabel('Quadrature');
    title(['Constelacao ',num2str(M(ii)),'-QAM']);
    hold on;
    ax = max(real(const_QAM))+sqrt(2)/2;
    line([0,0],[-ax,ax], 'Color', black);
    line([-ax,ax],[0,0],'Color', black);
    % saveas(h,['fig/',num2str(M(ii)),'_QAM_plot'],'pdf');
end 

pause;
close all;

const_4QAM = const_MQAM(M(1),distancia_QAM(1));
const_16QAM = const_MQAM(M(2),distancia_QAM(2));
const_64QAM = const_MQAM(M(3),distancia_QAM(3));


%% Grey Code

% Em_No = 30;     % - Razao sinal-ruido

% for jj = 1:M(ii)
%     bin_alfabeto = de2bi(0:jj-1);
%     text(real(const_QAM(jj)),imag(const_QAM(jj))+0.3, num2str(bin_alfabeto(jj,:)));
% end