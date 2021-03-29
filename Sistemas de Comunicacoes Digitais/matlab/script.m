%% Trabalho de SCD 
% Testes
% 2021/03/24 - Lucas Abdalah
%

close all; clearvars; clc; % Clear the matlab ambient

% addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema1'
addpath 'problema1' % Local path
addpath 'problema4' % Local path
addpath 'problema4\parte1' % Local path

M = 4;         % - Numero de simbolos da constelacao
% L = 4096;       % - Tamanho da sequencia (bits)
K = log2(M);    % - Numero de bits/simbolo
% N = L/K;        % - Numero de simbolos a serem transmitidos

% s_m = randi([0 1],L,1)';   % - Mensagem a transmitir

% Generating a gray const.
% G = gray_const(M,true);

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


for ii=1:sqrt(M)
    phi_i(1,ii) = ((2*ii - 1)*pi)/M(1);
end

const_PSK = zeros(1,M); % - Vetor que recebe o alfatbeto QAM.

k = 1;
for ii = 1:sqrt(M)      % bi - moving on the imaginary axis  
        const_PSK(1,ii) = cos(phi_i(1,jj)) + 1j*sin(phi_i(1,ii));
end

%% Rebater os simbolos
const_PSK = sqrt(E_s)*const_PSK;

pause
% This part is to display the symbol's bits (in Gray code)
% Organize the M-PSK modulation to sort the binary constellation and convert it into Gray encoding;
mat_PSK = reshape(const_PSK,[sqrt(M),sqrt(M)])';
mat_PSK(2:2:end,:) = fliplr(mat_PSK(2:2:end,:));
aux = mat_PSK.';
const_PSK = flip(aux(:)');
