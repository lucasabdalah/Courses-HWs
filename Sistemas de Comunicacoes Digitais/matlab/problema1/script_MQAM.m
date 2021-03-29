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
plot_figures = true;    % - Plotr ou nao
save_figures = false;   % - Salvar o plot ou nao
teste_map = false;      % - Entrar ou nao no teste das funcoes map e demap

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
        h(ii) = scatter_MQAM(const_QAM,true,false,true);    
        if save_figures == true
            saveas(h(ii),['fig/',num2str(M(ii)),'_QAM_plot'],'pdf');
        end
    end 
end 
disp('---------------------------------------')
pause
clc;
% close all;

%% Creating the constellation for each M-QAM
const_4QAM = const_MQAM(M(1),d);
const_16QAM = const_MQAM(M(2),d);
const_64QAM = const_MQAM(M(3),d);

%% Testar mapeador e demapeador com constelacao 4-QAM
if teste_map == true
    K = log2(M(1));    % - Numero de bits/simbolo
    L = 100*M(1)*K;       % - Tamanho da sequencia (bits)
    N = L/K;        % - Numero de simbolos a serem transmitidos

    %% cadeia de bits enviados
    s_m = randi([0 1],L,1)';    % - Mensagem a transmitir
    %% cadeia de bits a serem recebidos
    y_m = zeros(1,L);           % - Vetor vazio para receber a mensagem
    %% Simbolos a serem mapeados e demapeados
    symb_map = zeros(1,N);      % - Vetor vazio para receber os simbolos mapeados
    symb_demap = zeros(1,N);    % - Vetor vazio para receber os simbolos demapeados
    
    n=0;    % - Index for symbols iteration
    for ii = 1:K:L
        n=n+1;   
        bits = s_m(ii:ii+K-1);
        symb_map(n) = mapping_MQAM(bits);
        [symb_demap(n),y_m(ii:ii+K-1)]=demapping_MQAM(symb_map(n),M,d);
    end 
    disp(['message length: ', num2str(L)]);
    disp(['bit error: ', num2str(sum(xor(s_m,y_m)))]);
end

%% Exemplo com constelacao 16-QAM
if teste_map == false;
    example = 2;
    Nsymbols = 100;
    figure;
    scatter_MQAM(const_MQAM(M(example),d),true,false,true);
    hold on 

    for ii = 1:Nsymbols
        %% Signal coef
        a = randn(1)/sqrt(1);
        b = randn(1)/sqrt(1);
        r = a + 1j*b;
        %% Demapping 
        [symb,bits]=demapping_MQAM(r,M(example),d);
        %% Mapping the obtained bits;
        symb = mapping_MQAM(bits);
        %% Scatter Plot
        scatter_MQAM(r,false,true,false); 
        %% Linha da distancia Euclidiana entre o simbolo da constelacao e o demapeado
        line([real(symb),real(r)],[imag(symb),imag(r)],...
        'Color','black',...
        'LineStyle','--');
    end
    hold off

end;