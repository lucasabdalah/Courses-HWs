%% Trabalho de SCD 
% ----------------------------
% Problema 3 - Modulacao M-PSK
% script_MPSK.m
% 2021/03/26 - Lucas Abdalah
%

close all; clearvars; clc; % Clear the matlab ambient

addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab'
%% General Setup
plot_figures = false;    % - Plotr ou nao
save_figures = false;   % - Salvar o plot ou nao
teste_map = true;      % - Entrar ou nao no teste das funcoes map e demap

%% Problema 4 - Calculo da energia media (E_m) e da distancia entre os simbolos
M = [4, 8];    % - Numero de simbolos da constelacao
E_g = 1;            % - Energia do pulso de transmissao

energia_media_PSK = zeros(1,size(M,2));
distancia_PSK = zeros(1,size(M,2));

for ii = 1:size(M,2)
    energia_media_PSK(1,ii) = energia_MPSK(M(ii),E_g);
    distancia_PSK(1,ii) = d_MPSK(M(ii),energia_media_PSK(1,ii));
    fprintf('--------------------- Modulacao %1d-PSK --------------------- \n', M(ii));
    fprintf('Energia (E_media) = %1.2d \n', energia_media_PSK(1,ii));
    fprintf('Energia de bit (E_media_bit) = %1.2d \n', (energia_media_PSK(1,ii))/(3*log2(M(ii))));    
    fprintf('Distancia (d) = %1.2d \n\n', distancia_PSK(1,ii));
end
fprintf('------------------------------------------------------------\n');
% Since the distances are equal, let's define a standard for all cases;
d=distancia_PSK;
E_s = energia_media_PSK;

%% Plotando Constelacoes
if plot_figures == true
    for ii = 1:size(M,2)
        fprintf('>>Pressione qualquer tecla para ver a constelacao %1d-PSK \n', M(ii));
        pause;
        const_PSK = const_MPSK(M(ii),E_s(ii));
        figure;
        h(ii) = scatter_MPSK(const_PSK,true,false,true);
        if save_figures == true
            saveas(h(ii),['fig/',num2str(M(ii)),'_PSK_plot'],'pdf');
        end
    end 
end 
disp('---------------------------------------')
pause
clc;

%% Creating the constellation for each M-PSK
const_4PSK = const_MPSK(M(1),E_s(1));
const_8PSK = const_MPSK(M(2),E_s(2));

%% Testar mapeador e demapeador com constelacao 4-QPSK
if teste_map == true
    example = 2;
    K = log2(M(example));    % - Numero de bits/simbolo
    L = 100*M(example)*K;       % - Tamanho da sequencia (bits)
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
        symb_map(n) = mapping_MPSK(bits);
        [symb_demap(n),y_m(ii:ii+K-1)]=demapping_MPSK(symb_map(n),M(example),E_s(example));
    end 
    disp(['message length: ', num2str(L)]);
    disp(['bit error: ', num2str(sum(xor(s_m,y_m)))]);
end

%% Exemplo com constelacao 8-PSK
if teste_map == false;
    example = 2;
    Nsymbols = 20;
    figure;
    scatter_MPSK(const_MPSK(M(example),E_s(example)),true,false,true);
    hold on 

    for ii = 1:Nsymbols
        %% Signal coef
        if ii < Nsymbols/2
            r = 3*E_s(example)*exp(1j*pi*rand(1))/2;
        else 
            r = -3*E_s(example)*exp(1j*pi*rand(1))/2;
        end
        %% Demapping 
        [symb,bits]=demapping_MPSK(r,M(example),E_s(example));
        %% Mapping the obtained bits;
        symb = mapping_MPSK(bits);
        %% Scatter Plot
        scatter_MPSK(r,false,true,false); 
        %% Linha da distancia Euclidiana entre o simbolo da constelacao e o demapeado
        line([real(symb),real(r)],[imag(symb),imag(r)],...
        'Color','black',...
        'LineStyle','--');
    end
    hold off

end;