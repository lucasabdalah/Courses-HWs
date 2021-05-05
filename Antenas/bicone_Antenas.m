%% Trabalho de Antenas
% 2021/04/10 - Lucas Abdalah
%
%
%% Limpar o ambiente do matlab
close all; clearvars; clc;
%% Parte 2 - Antena Biconica
%% General def
plotFrequency = 2000000000; % Define plot frequency 
freqRange = (1800:20:2200) * 1e6; % Define frequency range 
% Define antenna 
load('bicone015_antennaDesigner.mat')
antennaObject = bicone015_antennaDesigner;
% 3-d pattern
figure; pattern(antennaObject, plotFrequency) 
title('Length = 0.015')
pause
clear antennaObject; 
% close all; clc;
load('bicone15_antennaDesigner.mat')
% Define antenna 
antennaObject = bicone15_antennaDesigner;
% 3-d pattern
figure; pattern(antennaObject, plotFrequency)
title('Length = 0.15');