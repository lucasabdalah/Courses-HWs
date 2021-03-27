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
        h(ii) = scatter_MQAM(const_QAM,true);    
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
% close all;
pause

%% Creating the constellation for each M-QAM
const_4QAM = const_MQAM(M(1),d);
const_16QAM = const_MQAM(M(2),d);
const_64QAM = const_MQAM(M(3),d);
disp('---------------------------------------')
clc;


%% Testar mapeador e demapeador
M = 4;         % - Numero de simbolos da constelacao
K = log2(M);    % - Numero de bits/simbolo
L = 100*M*K;       % - Tamanho da sequencia (bits)
N = L/K;        % - Numero de simbolos a serem transmitidos

s_m = randi([0 1],L,1)';   % - Mensagem a transmitir
y_m = zeros(1,L);
symb_map = zeros(1,N);
symb_demap = zeros(1,N);
n=0;
for ii = 1:K:L
    n=n+1;   
    bits = s_m(ii:ii+K-1);
    symb_map(n) = mapping_MQAM(bits);
    [symb_demap(n),y_m(ii:ii+K-1)]=demapping_MQAM(symb_map(n),M,d);
end 
disp(['message length: ', num2str(L)]);
disp(['bit error: ', num2str(sum(xor(s_m,y_m)))]);


% const_16QAM

% example = 2;
% figure;
% scatter_MQAM(const_16QAM,true);
% hold on 

% for jj = 1:M(example)
%     gray_alfabeto = gray_const(M(example),false);
%     str = {strjoin(string(gray_alfabeto(jj,:)))};
%     text(real(const_16QAM(jj)),imag(const_16QAM(jj))+(example/6),...
%     str,...
%     'FontSize', 6,...
%     'HorizontalAlignment', 'center');
% end
% % Draw a circle radius = d
% % center = [real(const_4QAM(jj)) imag(const_4QAM(jj))];
% % radii = sqrt(2)/2;
% % viscircles(center,radii,'Color','b');
% % axis square

% for ii = 1:4
%     % M = 16;
%     a = randn(1)/sqrt(1);
%     b = randn(1)/sqrt(1);
%     r = a + 1j*b;
%     [symb,bits]=demapping_MQAM(r,M(example),d);
%     disp(['ii=' num2str(ii)]);
%     bits
%     disp('---------------------------------------')
%     symb = mapping_MQAM(bits);
%     scatter(real(r),imag(r), 'filled',...
%     'MarkerEdgeColor',[1,0,0],... 
%     'MarkerFaceColor',[1,0,0]);
%     text(real(r),imag(r)-(example/15), num2str(ii),...
%     'FontSize', 6,...
%     'HorizontalAlignment', 'center');

%     line([real(symb),real(r)],[imag(symb),imag(r)]);
% end

% hold off

% bits = [1,1];
%% --------------------------------------------------------------------
pause;

% Em_No = 30;     % - Razao sinal-ruido