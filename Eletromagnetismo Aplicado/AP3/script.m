%%%%%%%%% METODO DOS MOMENTOS %%%%%%%%%
%%%%% Parte 1
% Definicao da Geometria, Aplicacao da Discretizacao e Calculo da Integral
%%%%% Parte 2
% Definicao da Geometria, Aplicacao da Discretizacao e Calculo da Integral
% Lucas de S. Abdalah - Eletromagnetismo Aplicado

% Inicializacao do Programa
clc; clear; close all;

%%%% Dados para Resolucao da questao
V0=2; % Volts
a=8e-2; % metros
b=7e-2; % metros
d=5e-3; % metros
%%%% Dados para Resolucao da questao
e0=8.8541e-12; % Constante epsilon zero

%%%%%%%%%%%%%%%%%%%% Parte 1 %%%%%%%%%%%%%%%%%%%%

% Discretizacao do Sistema
N=2.5e3; N_x2=2*N; 
M=sqrt(N); 
dx=a/M;  dy=b/M;  dL=dx;   
% Discretizacao do Sistema


%%% Adaptacao do Problema 15 da sessao Metodo dos Momentos do Sadiku vol. 3

%Cálculos dos coeficientes
aux=0; 

for p1_aux=1:2  

    for p2_aux=1:M  

        for p3_aux = 1:M 
            aux=aux+1; 
    
            X(aux)=dx*(p2_aux-0.5); 
    
            Y(aux)=dy*(p3_aux-0.5); 
        end  
    end  
end  
 
for p1_aux=1:N 
    
    Z(p1_aux)=0.0; 
    
    Z(p1_aux+N)=d; 

end 

for i=1:N_x2 
    
    for j=1:N_x2 
    
        if (i==j) 
    
            A(i,j) = dL*0.8814/(pi*e0); 
    
        else 
    
            R= sqrt((X(i)-X(j))^2 + (Y(i)-Y(j))^2+(Z(i)-Z(j))^2); 
    
            A(i,j) = dL^2/(4.*pi*e0*R); 
    
        end 
    end 
end 

% Matriz B descrita na introducao do trabalho 
for aux=1:N 

    B(aux)=2; 

    B(aux+N)=-4; 

end 
 
RHO=inv(A)*B'; % Resolvendo o sistema descrito

soma=0; 
 
n_local=1000; 
 
for i=1:n_local
    
    soma=soma+RHO(i);  

end 

% Caracteristcas da Malha
X=linspace(-4,4,M); Y=linspace(-7,7,M); Z=zeros(M);
                       

% Calculo de carga na superficie 
for i =1:M 
    
    for j=1:M 
    
        Z(i,j)=dx*dy*RHO(j+(i-1)*M); 
    end    
end 
 
Z_inferior=zeros(M); 
 
for i=M+1:M+M 
    
    for j=1:M 
    
        Z_inferior(i-M,j)=dx*dy*RHO(j+(i-1)*M); 
    
    end 
end 
%%% Adaptacao do Problema 15 da sessao Metodo dos Momentos do Sadiku vol. 3


% Plot das superficies no espaço
surf(X,Y,Z); 
hold on 
surf(X,Y,Z_inferior); 
title({'Problema 1'; 'Geometria Aproximada do Problema'}); % Titulo do Grafico
xlabel('x'); % Titulo do eixo x 
ylabel('y'); % Titulo do eixo y
zlabel('z'); % Titulo do eixo z


% COMENTARIO
hold off
pause
disp('Pressione qualquer tecla para seguir (ver os graficos).');
clc; close all;
% COMENTARIO



%%%%%%%%%%%%%%%%%%%% Parte 2 %%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

%%%% Dados para Resolucao da questao
v1=2; % Volts
v2=-5; % Volts
v0=v1-v2; % Volts
%%%% Dados para Resolucao da questao

M=40;
I= 5000;

% Caracteristcas da Malha
X=linspace(-6,6,2*M); Y=linspace(-4,4,2*M); 

% Matriz a ser preenchida por potenciais
V1=zeros(2*M);

% Matriz preenchida pelo metodo iterativo para calcular o potencial desconhecido
for i=1:2*M;
    
    for j=1:2*M;
    
        if (i>=5 && i<=75 && j>= 15 && j <= 65);
    
            V1(i,j)=v0;
    
        end
    end
end

for k=1:I 
    
    for i=2:2*M-1
    
        for j=2:2*M-1
    
            if i>= 2 && V1(i,j)~=v0
    
                V1(i,j) = (V1(i-1,j) + V1(i+1,j) + V1(i,j-1) + V1(i,j+1))/4;
    
            end
        end
    end
end

[Ex, Ey] = gradient(V1); % Obtem campo vetorial a partir do gradiente do Potencial calculado
figure(2)
contour(X,Y,V1,10) % Equipotenciais
hold on
quiver(X,Y,-Ex./sqrt(Ex.^2 + Ey.^2),-Ey./sqrt(Ex.^2 + Ey.^2)) % Campo Eletrico normalizado

title({'Problema 1'; 'Campo Eletrico e Linhas Equipotenciais'}); % Titulo do Grafico
xlabel('x (cm)'); % Titulo do eixo x 
ylabel('y (cm)'); % Titulo do eixo y
disp('Aperte enter para fechar tudo.');

% COMENTARIO
hold off
pause
clc; close all;
% COMENTARIO

% Caracteristcas da Malha
Z=linspace(-3.25,3.25,2*M); Y=linspace(-6,6,2*M); 

% Matriz a ser preenchida por potenciais
V1=zeros(2*M);

% Aproximacao de valores discretos nas coordenadas utilizadas para o problema
x1= round(2/6*80); x2=round(4/6*80);
y1= round(2/12*80); y2= round(10/12*80);


for i=1:2*M;
    
    for j=1:2*M;
    
        if j>=y1 && j<=y2;
    
            if i== x1;
    
                V1(i,j)=v1;
    
            end
    
            if i== x2;
    
            V1(i,j)=v2;
    
            end
        end
    end
end

for k=1:I 
    
    for i=2:2*M-1
    
        for j=2:2*M-1
    
            if i>= 2 && V1(i,j)~=v1 && V1(i,j)~=v2
    
                V1(i,j) = (V1(i-1,j) + V1(i+1,j) + V1(i,j-1) + V1(i,j+1))/4;
    
            end
        end
    end
end

figure(1)
contourf(Y,Z,V1,30)%,'edgecolor','none')%determina as linhas equipotenciais associadas ao campo escalar e plota
colorbar %fixa no plot a barra de cores que mapeia os valores das curvas potenciais

title({'Problema 1'; 'Distribuicao de Potencial'}); % Titulo do Grafico
xlabel('y'); % Titulo do eixo x 
ylabel('z'); % Titulo do eixo y
set(gca,'YDir','reverse')
disp('Aperte enter para fechar tudo.');

% COMENTARIO
hold off
pause
clc; close all;
% COMENTARIO

[Ey, Ez] = gradient(V1); %calcula as componentes do campo vetorial
figure(1)
contour(Y,Z,V1,10) %determina as linhas equipotÊnciais associadas ao campo escalar                 
hold on
quiver(Y,Z,-Ey./sqrt(Ey.^2 + Ez.^2),-Ez./sqrt(Ey.^2 + Ez.^2)) %normaliza o campo elétrico e então o plota
set(gca,'YDir','reverse')

title({'Problema 1'; 'Campo Eletrico e Linhas Equipotenciais'}); % Titulo do Grafico
xlabel('y (cm)'); % Titulo do eixo x 
ylabel('z (cm)'); % Titulo do eixo y
disp('Aperte enter para fechar tudo.');

