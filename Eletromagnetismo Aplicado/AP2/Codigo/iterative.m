%%%%%%%%% METODO ITERATIVO %%%%%%%%%
%%%%% Parte 1
% Script para calculo da distribuição de potencial utilizando o método das diferenças finitas iterativo.
%%%%% Parte 2
% Plota Campo Vetorial do campo elétrico associado à distribuição de potencial calculada. 
% Lucas de S. Abdalah - Eletromagnetismo Aplicado

% Inicializacao do Programa
clc; clear; close all;

% Dados para resolucao da questaoV0=11; % Volts
x_1=3; % cm
x_2=5; % cm 
x_3=7; % cm
x_4=9; % cm
y_1=2; % cm
y_2=0; % cm
d=2; % cm

% Linha de Potencial - Sequencia dos pontos	
% (x_1,y_1); (x_2,y_1); (x_2,y_2); (x_3,y_2); (x_3,y_1); (x_4,y_1)

%%%%% Parte 1
N=50; % Quantidade de pontos do intervalo
h=10; % cm % Limite da Geometria
X=linspace(0,h,N); % Geometria do eixo X
Y=linspace(0,h,N); % Geometria do eixo Y

% Adaptacao das coordenadas para o caso de N
x1=floor((x_1*N)/h + 1);
x2=floor((x_2*N)/h + 1);
x3=floor((x_3*N)/h + 1);
x4=floor((x_4*N)/h + 1);
y1=floor((y_1*N)/h + 1);
y2=floor((y_2*N)/h + 1);

% Potencial e Tolerancia
V=zeros(N); % Matriz a ser preenchida por potenciais
tol=ones(N)*1e-6; % Vetor de tolerancia para criterio de parada

maxtol=1e6; % Maximo de Iteracoes

% Metodo Iterativo
for aux=1:maxtol

	verifica=V;
	
	for i=2:N-1
		
		for j=2:N-1
			
			if (i==y2+1 && j>=x2 && j<=x3) || (i==y1 && j>=x1 && j<=x2) || (i<=y1 && i>=y2 && j==x2) ||  (i<=y1 && i>=y2 && j==x3) || (i==y1 && j>=x3 && j<=x4)
			
				V(i,j)=V0;
				continue;
				
			else
				V(i,j) = (V(i-1,j) + V(i+1,j) + V(i,j-1) + V(i,j+1))/4;	
			end
		end
	end
	
	if V-verifica<tol

		disp('Numero de iteracoes necessarias:');
		disp(aux);
		break;

	else
		continue;

	end
end

disp('Pressione qualquer tecla para ver os graficos.');

pause;

%%%%% Parte 2
X=X(2:end);Y=Y(2:end);V=V(2:end,2:end); % Ajuste para plot mais fidedigno a distribuicao de potencial

figure(1)
contour(X,Y,V,1000)  % Distribuicao de potencial associada ao campo escalar
colorbar % Introduz a barra de cores pra mapear as curvas equipotencias. 
title({'Problema 1'; 'Distribuicao de Potencial'}); % Titulo do Grafico
xlabel('x (cm)'); % Titulo do eixo x 
ylabel('y (cm)'); % Titulo do eixo y

[Ex, Ey] = gradient(V); % Calculo do Gradiente, associado ao campo escalar

figure(2)
contour(X,Y,V,10) % Determina as linhas equipotenciais associadas ao campo escalar e plota
hold on
quiver(X,Y,-Ex,-Ey) % Plota os vetores (EX,EY) nos pontos (X,Y) no espaço bidimensional
colorbar
title({'Problema 1'; 'Gradiente e Linhas Equipotenciais'}); % Titulo do Grafico
xlabel('x (cm)'); % Titulo do eixo x 
ylabel('y (cm)'); % Titulo do eixo y

figure(3)
contour(X,Y,V,10) % Determina as linhas equipotenciais associadas ao campo escalar e plota
hold on
quiver(X,Y,-Ex./sqrt(Ex.^2 + Ey.^2),-Ey./sqrt(Ex.^2 + Ey.^2)) % plota os vetores (EX,EY) nos pontos (X,Y) no espaço bidimensional 
colorbar
title({'Problema 1'; 'Gradiente e Linhas Equipotenciais (Normalizado)'}); % Titulo do Grafico
xlabel('x (cm)'); % Titulo do eixo x 
ylabel('y (cm)'); % Titulo do eixo y

% Finalizacao do Programa
hold off

disp('Aperte enter para fechar tudo.');
pause;

clc; close all; clear all