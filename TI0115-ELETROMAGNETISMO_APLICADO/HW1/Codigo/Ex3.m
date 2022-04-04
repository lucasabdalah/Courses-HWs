% 	Plotar o campo vetorial

clc; % limpa a tela

[X,Y] = meshgrid(-4:.4:4); % gera a grade de pontos - no caso bidimensional 

a = 4.*Y; b = 1.8.*X; c = 8;

EX = a./(X.^2 + Y.^2 + c); % componente Ex do campo vetorial E
EY = b./(X.^2 + Y.^2 + c); % componente Ey do campo vetorial E
hold on
quiver(X,Y,EX,EY)
colormap hsv

hold off
title({'Problema 3'; 'Campo Eletrico'}); % Titulo do Grafico
xlabel('Eixo X') % Titulo do eixo x 
ylabel('Eixo Y') % Titulo do eixo y 