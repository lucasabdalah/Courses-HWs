%% AP2 de PES
% Questao 2. c)
%
% item_c.m
%
% 2021/08/23 - Lucas Abdalah

close all; clearvars; clc; % Clear the matlab ambient

%% General Setup
pi_0 = 0.3; % Prob a priori P(H_0) 
pi_1 = 0.7; % Prob a priori P(H_1) 
C_10 = 1; % Custo de Decisao 1 para H_0 correta
C_01 = 2; % Custo de Decisao 0 para H_1 correta

%% Valores de a
a = 0:1e-2:1; % a variando de 0 a 1 

limiar = (2*log(a*C_10./((1-a)*C_01)) + 1)./2;

%% Funcoes Risco
R_0 = C_10*(qfunc(limiar));
R_1 = C_01*(qfunc(1-limiar));
R = 2*a.*R_0 + (1-a).*R_1;

%% Plot 
h = figure();
plot(a, R,...
'Marker','x',...
'LineWidth', 1.5,...
'LineStyle', '-');
title('Minimax: Risco');
xlabel('Probabilidade, $a$','interpreter','latex');
ylabel('Risco, $R$','interpreter','latex');
grid on
xticks(0:.1:1);
[M_R I_R] = max(R);
x_arrow = [0.53,0.53];
y_arrow = [0.8, 0.9];
a_arrow = annotation('textarrow',x_arrow,y_arrow,'String',['a = ',num2str(a(I_R))])
% saveas(h,'item_c_Minimax.png');