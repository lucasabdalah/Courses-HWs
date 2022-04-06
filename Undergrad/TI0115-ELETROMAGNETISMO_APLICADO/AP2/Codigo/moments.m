%%%%%%%%% METODO DOS MOMENTOS %%%%%%%%%
%%%%% Parte 1
% Calculo das matrizes iniciais A e B para método dos momentos
%%%%% Parte 2
% Calculo do Potencial V e plot campo elétrico associado à distribuição de potencial calculada. 
% Lucas de S. Abdalah - Eletromagnetismo Aplicado

% Inicializacao do Programa
clc; clear; close all;

% Dados para resolucao da questao
eps0=8.8541e-12; %permissividade no vácuo. 
V0=11; % Volts
a=3e-6; 
b=3.85;
L=5e-2; 
N=100; % Discretizacao do filamento 
D=L/N; % Gradacao44              

% Criacao da matriz A detalhada no trabalho
A = zeros(N); % Matriz que recebe a geometria
I = 1:N;
Y = D*(I-0.5);

%%%%% Parte 1
for i = 1:N
 
    for j = 1:N
 
        if(i~=j)
 
            A(i,j) = D/abs(Y(i)-Y(j));
 
        else
 
            A(i,j) = 2*log(D/a);
 
        end
    end
end

% Criacao da Matriz B com os dados para resolucao
B = b*4*pi*eps0*V0*ones(N,1);

Rho = inv(A)*B; % Rho - Distriubicao de potencial

Qd = Rho*D; %Carga da discretizacao

% Dados fornecidos para resolucao da questao
k = 8.987E9; % Constante
x1 = 3e-2; % m
y1 = 2e-2; % m 
d = 2e-2; % m
disc = N/5; % Discretizacao

%%%%% Parte 2
% Matriz $\rho$ descrita no trabalho
p = [linspace(x1,x1+d,disc),... % eixo X
    linspace(x1+d,x1+d,disc),...
    linspace(x1+d,x1+2*d,disc),...
    linspace(x1+2*d,x1+2*d,disc),...
    linspace(x1+2*d,x1+3*d,disc);...
    
    linspace(y1,y1,disc),...    % eixo Y
    linspace(y1,y1-d,disc),...
    linspace(y1-d,y1-d,disc),...
    linspace(y1-d,y1,disc),...
    linspace(y1,y1,disc)].';
Q = Qd; % Elementos de Carga

% Grid para calculo do potencial
[X,Y] = meshgrid(linspace(0,.1,100),linspace(0,.1,100));  % Metros

% Matriz dos Potenciais a ser preenchida
V = zeros(size(X));  

% Preenche o potencial de cada carga
for aux = 1:numel(Q) 
    V = V + k * Q(aux) ./ sqrt(abs(p(aux,1)-X).^2+abs(p(aux,2)-Y).^2);
end

disp('Pressione qualquer tecla para ver os graficos.');

pause;

figure(1)
surf(peaks)
hContour = contourf(X,Y,V,'LevelList',[0 quantile(V(:),50)],'edgecolor','none'); % Distribuicao de potencial associada ao campo escalar
hColorbar = colorbar;  % Introduz a barra de cores pra mapear as curvas equipotencias. 
ylabel(hColorbar,'Potencial elétrico (V)')
title({'Problema 2'; 'Distribuicao de Potencial'}); % Titulo do Grafico
xlabel('x (m)'); % Titulo do eixo x 
ylabel('y (m)'); % Titulo do eixo y

% Calculo do Gradiente da Funcao
[Ex,Ey] = gradient(-V);

validColumns = all(isfinite(Ex) & isfinite(Ey)); % Atribui apenas colunas que 'convergem'

figure(2)
hContour = contourf(X,Y,V,'LevelList',[0 quantile(V(:),20)],'edgecolor','none'); % plota os vetores (EX,EY) nos pontos (X,Y) no espaço bidimensional 
hold on
colorbar
hLines = streamslice(X(:,validColumns),Y(:,validColumns),Ex(:,validColumns),Ey(:,validColumns));
set(hLines,'Color','w'); % Define a cor branca pras linhas de campo
title({'Problema 2'; 'Gradiente e Linhas Equipotenciais'}); % Titulo do Grafico
xlabel('x (m)'); % Titulo do eixo x 
ylabel('y (m)'); % Titulo do eixo y

% Finalizacao do Programa
hold off

disp('Aperte enter para fechar tudo.');
pause;

clc; close all; clear all