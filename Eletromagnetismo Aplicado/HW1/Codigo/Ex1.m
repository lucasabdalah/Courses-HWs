% 	Plotar as superfícies equipotenciais e o campo vetorial
%   associados ao campo escalar

clc;  % limpa a tela

[X,Y] = meshgrid(-1:.1:1);  % gera a grade de pontos - no caso bidimensional 

a = 2.0; b = 2.2; c = 3.*X; d = 4.*Y; % constantes fornecidas pelo problema

Ve = (c+d).*exp(a.*X.^2 + b.*Y.^2);	% função que define o campo escalar, usa a grade definida pela meshgrid

[DX,DY] = gradient(Ve,.1,.1); % calcula as componentes do campo vetorial associado ao campo escalar 
contour(X,Y,Ve) % determina as linhas equipotenciais associadas ca campo escalar e plota
hold on

DX = -1.*DX; % troca o sinal para calcular E - - grad
DY = -1.*DY; % troca o sinal para calcular E - - grad

quiver(X,Y,DX,DY) % plota os vetores (DX,DY) nos pontos (X,Y) no espaço bidimensional 
colormap hsv
hold off

title({'Problema 1'; 'Gradiente e Linhas Equipotenciais'}); % Titulo do Grafico
xlabel('Distancia (m)'); % Titulo do eixo x 
ylabel('Campo Eletrico (E)'); % Titulo do eixo y

figure	
[X,Y] = meshgrid(-1:.1:1); % gera a grade de pontos - no caso bidimensional 
contour(X,Y,Ve) % determina as linhas equipotenciais associadas ca campo escalar e plota
hold on

DX = DX./sqrt(DX.^2 + DY.^2); % componente Ex do campo vetorial E
DY = DY./sqrt(DX.^2 + DY.^2); % componente Ey do campo vetorial E

quiver(X,Y,DX,DY) % plota os vetores (DX,DY) nos pontos (X,Y) no espaço bidimensional 
colormap hsv
hold off

title({'Problema 1'; 'Gradiente e Linhas Equipotenciais'}); % Titulo do Grafico
xlabel('Distancia (m)'); % Titulo do eixo x 
ylabel('Campo Eletrico (E)'); % Titulo do eixo
grid on; % Mostra as grades