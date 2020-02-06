%%%%%%%%% MATRIZ BANDA %%%%%%% Parte 1
% Script para calculo da distribuição de potencial utilizando o método da matriz de banda.
%%%%% Parte 2
% Plota Campo Vetorial do campo elétrico associado à distribuição de potencial calculada. 

N=5;

% Matriz de banda usandi a referida aproximacao referente geomentria

A = [-4, 1,  0,  1,	zeros(1,N);	% 1
     1, -4, 1,  0, 1,  zeros(1,4);	% 2
     0,  1, -4,  0, 0,  1,  0,  0,  0;	% 3
     1,  0,  0, -4, 1,  0,  1,  0,  0;	% 4
     0,  1,  0,  1, -4,  1,  0,  1,  0;	% 5
     0,  0,  1,  0,  1, -4,  0,  0,  1;	% 6
     0,  0,  0,  1,  0,  0, -4,  1,  0;	% 7
     0,  0,  0,  0,  1,  0,  1, -4,  1;	% 8
     0,  0,  0,  0,  0,  1,  0,  1, -4]	% 9

% Lado direito da igualdade
B = [0;
	-16;
	-4;
	0;
	0;
	0;
	0;
	-4;
	-16];

A_i = inv(A); %Calculo da Inversa

V = A_i*B; % Obtendo a solucao V do sistema linear

% Reconstrucao da matriz original, com as fronteiras da geometria
v = [zeros(1,N);
     0, V(1), V(2), V(3), 0;
     0, V(4), V(5), V(6), 0;
     0, V(7), V(8), V(9), 0;
     			zeros(1,N)];


disp('Pressione qualquer tecla para ver os graficos.');

pause
     
%%%%% Parte 2
X = linspace(0,10,N); Y = linspace(0,8,N); % Cria os pontos da geometria

figure(1)
contour(X,Y,v,200) % Distribuicao de potencial associada ao campo escalar
colorbar % Introduz a barra de cores pra mapear as curvas equipotencias.
set(gca,'YDir','reverse')

title({'Problema 2'; 'Distribuicao de Potencial'}); % Titulo do Grafico
xlabel('x (cm)'); % Titulo do eixo x 
ylabel('y (cm)'); % Titulo do eixo y

[Ex, Ey] = gradient(v);  % calcula as componentes do campo vetorial associado ao campo escalar.

figure(2)
contour(X,Y,v,10) % Distribuicao de potencial associada ao campo escalar
hold on
quiver(X,Y,-Ex./sqrt(Ex.^2 + Ey.^2),-Ey./sqrt(Ex.^2 + Ey.^2)) % plota os vetores (EX,EY) nos pontos (X,Y) no espaço bidimensional
set(gca,'YDir','reverse')

title({'Problema 2'; 'Gradiente e Linhas Equipotenciais (Normalizado)'}); % Titulo do Grafico
xlabel('x (cm)'); % Titulo do eixo x 
ylabel('y (cm)'); % Titulo do eixo y


disp('Aperte enter para fechar tudo.')
pause

clc; close all; clear all