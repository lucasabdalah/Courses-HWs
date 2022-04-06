%% Trabalho de SCD
% 2021/03/24 - Lucas Abdalah
%

close all; clearvars; clc; % Clear the matlab ambient

% addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema1'
addpath 'problema1' % Local path
addpath 'problema2' % Local path
addpath 'problema3' % Local path

script_MQAM;
disp('>> pressione qualquer botao para continuar para questao 2');
pause;
close all; clearvars; clc; % Clear the matlab ambient

script_prob;
disp('>> pressione qualquer botao para continuar para questao 3 - parte 1');
pause;
close all; clearvars; clc; % Clear the matlab ambient

script_AWGN;
disp('>> pressione qualquer botao para continuar para questao 3 - parte 2');
pause;
close all; clearvars; clc; % Clear the matlab ambient

script_teoricaxAWGN;
disp('>> pressione qualquer botao para continuar para questao 4 - parte 1');
pause;
close all; clearvars; clc; % Clear the matlab ambient

addpath 'problema4\parte1' % Local path
addpath 'problema4\parte2' % Local path
addpath 'problema4\parte3' % Local path
rmpath('problema1') % Local path
rmpath('problema2') % Local path
rmpath('problema3') % Local path

script_MPSK;
disp('>> pressione qualquer botao para continuar para questao 4 - parte 2');
pause;
close all; clearvars; clc; % Clear the matlab ambient

script_prob;
disp('>> pressione qualquer botao para continuar para questao 4 - parte 3');
pause;
close all; clearvars; clc; % Clear the matlab ambient

addpath 'problema1' % Local path
addpath 'problema2' % Local path
addpath 'problema3' % Local path
addpath 'problema4\parte1' % Local path
addpath 'problema4\parte2' % Local path
addpath 'problema4\parte3' % Local path

script_teoricaxAWGN;
disp('>> pressione qualquer botao para continuar para questao 5');
pause;
close all; clearvars; clc; % Clear the matlab ambient


