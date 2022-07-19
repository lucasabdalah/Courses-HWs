%% Função Main
% [TIP7188 - Filtragem Adaptativa]
% Author: Lucas Abdalah
%
%
% main.m

clearvars;
close all;
clc; pause(0.1)

% publish('main.m', 'pdf');
% publish('filter_hw.m', 'pdf');

% filter_hw.hw2p5(false);

Rx = [1,0.5; 0.5,1];
Pxd = [0.5; 0.25];
invRx = inv(Rx);
w = invRx*Pxd;