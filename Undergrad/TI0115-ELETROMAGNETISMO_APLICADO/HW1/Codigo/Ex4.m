% 	Plotar o campo vetorial

clc; % limpa a tela

% X = 0; % Fornecido pelo problema

[Y,Z] = meshgrid(0:.5:3,-3:.5:3); % Gera a grade de pontos para transformacaoo de coordenadas

raio = sqrt(Y.^2 + Z.^2); % Conversao da coordenada raio
theta = acos(Z./raio); % Conversao da coordenada theta
phi = pi/2;	% Resultado obtido a partir do enunciado

teste = (4.*sin(theta));
teste1 = (3.5.*cos(theta));

Araio = teste./(raio.^2 + 6.5); % Referente ao campo na direcao do raio
Atheta = teste1./(raio.^2 + 6.5); % Referente ao campo na direcao de theta
Aphi = 0; % Referente ao campo na direcao phi

[X,Y,Z,Ax,Ay,Az] = VetEsfRet(raio,theta,phi,Araio,Atheta,Aphi); % Transformacao utilizando a rotina proposta

hold on
quiver(Y,Z,Ay,Az,'b') % plota os vetores (DX,DY) nos pontos (X,Y) no espaço bidimensional 
colormap hsv

[Y,Z] = meshgrid(-3:.5:0,-3:.5:3); % Gera a grade de pontos para transformacaoo de coordenadas

raio = sqrt(Y.^2 + Z.^2); % Conversao da coordenada raio
theta = acos(Z./raio); % Conversao da coordenada theta
phi = 3*pi/2;	% Resultado obtido a partir do enunciado

teste = (4.*sin(theta));
teste1 = (3.5.*cos(theta));

Araio = teste./(raio.^2 + 6.5); % Referente ao campo na direcao do raio
Atheta = teste1./(raio.^2 + 6.5); % Referente ao campo na direcao de theta
Aphi = 0; % Referente ao campo na direcao phi

[X,Y,Z,Ax,Ay,Az] = VetEsfRet(raio,theta,phi,Araio,Atheta,Aphi); % Transformacao utilizando a rotina proposta

quiver(Y,Z,Ay,Az,'b') % plota os vetores (DX,DY) nos pontos (X,Y) no espaço bidimensional 
colormap hsv

title({'Problema 4'; 'Campo Eletrico'}); % Titulo do Grafico
xlabel('Eixo Y') % Titulo do eixo x 
ylabel('Eixo Z') % Titulo do eixo y
grid on; % Mostra as grades

