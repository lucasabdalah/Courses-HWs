% 	Plotar o campo vetorial

clc; % limpa a tela

[ro,phi] = meshgrid(0.1:.05:0.6,0:.05:pi/6); % Gera a grade de pontos para transformacaoo de coordenadas

z = 0;	% Resultado obtido a partir do enunciado

Ve = (360.*phi)./pi; % Potencial definido pelo problema

% Gradiente calculado teoricamente (devido facilidade)
Gradraio = 0; % Referente ao campo na direcao do raio
Gradphi = -360./(ro.*pi); % Referente ao campo na direcao phi
Gradz = 0; % 

[X,Y,Z,Ax,Ay,Az] = VetCilRet(ro,phi,z,Gradraio,Gradphi,Gradz); % Transformacao utilizando a rotina proposta

hold on
quiver(X,Y,Ax,Ay) % plota os vetores (DX,DY) nos pontos (X,Y) no espaço bidimensional 
colormap hsv
[X,Y,Z] = CooCilRet(ro,phi,Ve);
contour(X,Y,Z) % linhas equipotenciais associadas ao campo escalar
hold off

title({'Problema 2'; 'Campo Eletrico'}); % Titulo do Grafico
xlabel('Eixo X') % Titulo do eixo x 
ylabel('Eixo Y') % Titulo do eixo y
grid on; % Mostra as grades