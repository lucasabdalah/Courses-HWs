%%%%%%%%% CALCULO DE POTENCIAL E CAMPO ELETRICO DE ESPIRA %%%%%%%%%
%%%%% Parte 1
% Definicao da Geometria, Aplicacao da Discretizacao e Calculo da Integral
%%%%% Parte 2
% Visualizacao da Geometria, Potencial e Campo Eletrico
% Lucas de S. Abdalah - Eletromagnetismo Aplicado

% Inicializacao do Programa
clc; clear; close all;

%%%%%%%%%%%%%%%%%%%% Parte 1 %%%%%%%%%%%%%%%%%%%%

	%%%% Dados para Resolucao da questao
	A=7e-2; % metro
	B=11e-2; % metro
	a=A/2; b=B/2;
	RHOL=5; % Densidade de Distribuicao de Carga proposta
	%%%% Dados para Resolucao da questao
	e0=8.854e-12; % Constante epsilon zero
	N=50; % Discretizacaoa da geometria

	% PLANO XY
	X=linspace(-a,a,N); Y=linspace(-b,b,N);

	% Definicao da discretizacao
	dx=a/N; dy=b/N;

	% Calculo da Integral
	x=a/2; y=b/2;
	somax=0; somay=0;

	for i = 0:N
		
		Dy = sqrt(x^2+(-a/2 + Y*dy).^2);
		
		somay = somay +1./Dy;
	end

	% Aplicando a equacao do Potencial
	Vy = (RHOL*somay*dy)/4*pi*e0;

	for i = 0:N

		Dx = sqrt(y^2+(-b/2 + X*dx).^2);

		somax = somax +1./Dx;
	end

	% Aplicando a equacao do Potencial
	Vx = (RHOL*somax*dx)/4*pi*e0;

	% Definido pela geometria
	Vc = 2*Vx +2*Vy;

	% Matriz de potenciais a ser preenchida
	V = zeros(N);

	% Aproximacao de valores discretos nas coordenadas utilizadas para o problema
	xc = round(((a/2)*N/a)+1); yc = round(((b/2)*N/b)+1);

	V(xc,yc) = 4;

	% Por diferenças finitas
	I = 1000;

	for k=1:I
		
		for i=2:N-1
		
			for j= 2:N-1
		
				V(i,j) = (V(i-1,j) + V(i+1,j) + V(i,j-1) + V(i,j+1))/4; %Calcula o potencial de um ponto que não foi pré-definido fazendo a média dos 4 pontos ao seu redor.
			end
		end
	end

%%%%%%%%%%%%%%%%%%%% Parte 2 %%%%%%%%%%%%%%%%%%%%

	%%%%% VISUALIZACAO GERAL
	% Vertices no eixo x
	X_espira=[a a -a -a a];
	% Vertices no eixo y
	Y_espira=[-b b b -b -b];
	% Vertices no eixo z
	Z_espira=zeros(1,5);

	figure1 = figure; % Create figure
	axes1 = axes('Parent',figure1); % Create axes
	hold(axes1,'on');
	fill3(X_espira,Y_espira,Z_espira,'y');
	grid on;
	title({'Geometria da Espira - 3D'});
	xlabel('x(m)'); 
	ylabel('y(m)'); 
	zlabel('z(m)');
	view(axes1,[-55.5 55.6]);

	% COMENTARIO
	hold off
	disp('Pressione qualquer tecla para ver os graficos.');
	pause
	clc; close all;
	% COMENTARIO

	%%%%%%%%%%%%% PLANO XY %%%%%%%%%%%%
	figure1 = figure; % Create figure
	axes1 = axes('Parent',figure1); % Create axes
	hold(axes1,'on');
	fill3(X_espira,Y_espira,Z_espira,'y');
	grid on;
	title({'Plano XY'});
	xlabel('x(m)'); 
	ylabel('y(m)'); 
	view(axes1,[0 90]);

	% COMENTARIO
	hold off
	disp('Pressione qualquer tecla para ver os graficos.');
	pause
	clc; close all;
	% COMENTARIO

	%%% VISUALIZACAO DO PLANO XY %%%
	figure
	contour(X,Y,V,8)
	grid on
	title({'Problema 1'; 'Linhas Equipotenciais'}); % Titulo do Grafico
	xlabel('x (m)'); % Titulo do eixo x 
	ylabel('y (m)'); % Titulo do eixo y

	% COMENTARIO
	hold off
	disp('Pressione qualquer tecla para ver os graficos.');
	pause
	clc; close all;
	% COMENTARIO

	contourf(X,Y,V,100,'edgecolor','none')%Determina as linhas equipotenciais associadas ao campo escalar e plota 
	colorbar %Fixa no plot a barra de cores que mapeia os valores das curvas potenciais 
	title({'Problema 1'; 'Distribuicao de Potencial no plano XY'}); % Titulo do Grafico
	xlabel('x (m)'); % Titulo do eixo x 
	ylabel('y (m)'); % Titulo do eixo y

	% COMENTARIO
	hold off
	disp('Pressione qualquer tecla para ver os graficos.');
	pause
	clc; close all;
	% COMENTARIO

	[Ex,Ey] = gradient(V);

	figure(2)
	contour(X,Y,V,8)
	hold on
	quiver(X,Y,-Ex,-Ey)

	title({'Problema 1'; 'Campo Eletrico (Nao Normalizado) e Linhas Equipotenciais'}); % Titulo do Grafico
	xlabel('x (m)'); % Titulo do eixo x 
	ylabel('y (m)'); % Titulo do eixo y
	disp('Aperte enter para fechar tudo.');

	% COMENTARIO
	hold off
	disp('Pressione qualquer tecla para ver os graficos.');
	pause
	clc; close all;
	% COMENTARIO

	%%% VISUALIZACAO DO PLANO XZ
	figure1 = figure; % Create figure
	axes1 = axes('Parent',figure1); % Create axes
	hold(axes1,'on');
	fill3(X_espira,Y_espira,Z_espira,'y');
	grid on;
	title({'Plano XZ'});
	xlabel('x(m)');
	zlabel('z(m)');
	view(axes1,[0 0]);


	% COMENTARIO
	hold off
	disp('Pressione qualquer tecla para ver os graficos.');
	pause
	clc; close all;
	% COMENTARIO

	figure1 = figure; % Create figure
	axes1 = axes('Parent',figure1); % Create axes
	contourf(X,Y,V,100,'edgecolor','none')%Determina as linhas equipotenciais associadas ao campo escalar e plota 
	colorbar %Fixa no plot a barra de cores que mapeia os valores das curvas potenciais 
	title({'Problema 1'; 'Distribuicao de Potencial no plano XZ'}); % Titulo do Grafico
	xlabel('x (m)'); % Titulo do eixo x 
	ylabel('y (m)'); % Titulo do eixo y
	zlabel('z(m)');
	view(axes1,[0 0.5]);

	% COMENTARIO
	hold off
	disp('Pressione qualquer tecla para ver os graficos.');
	pause
	clc; close all;
	% COMENTARIO

	%%% VISUALIZACAO DO PLANO YZ
	figure1 = figure; % Create figure
	axes1 = axes('Parent',figure1); % Create axes
	hold(axes1,'on');
	fill3(X_espira,Y_espira,Z_espira,'y');
	grid on;
	title({'Plano YZ'}); 
	ylabel('y(m)'); 
	zlabel('z(m)');
	view(axes1,[-90 0]);

	% COMENTARIO
	hold off
	disp('Pressione qualquer tecla para ver os graficos.');
	pause
	clc; close all;
	% COMENTARIO
	figure1 = figure; % Create figure
	axes1 = axes('Parent',figure1); % Create axes
	contourf(X,Y,V,100,'edgecolor','none')%Determina as linhas equipotenciais associadas ao campo escalar e plota 
	colorbar %Fixa no plot a barra de cores que mapeia os valores das curvas potenciais 
	title({'Problema 1'; 'Distribuicao de Potencial no plano YZ'}); % Titulo do Grafico
	xlabel('x (m)'); % Titulo do eixo x 
	ylabel('y (m)'); % Titulo do eixo y
	zlabel('z(m)');
	view(axes1,[90 0.5]);

% COMENTARIO
hold off
disp('Pressione para fechar tudo.');
pause
clc; close all; clear all;
% COMENTARIO