function Y = unifor_GenCauchy() %Função de geração da distribuição de Cauchy
alfa = 3; %Atribui valor 4 ao lambda
k = 10000; %10000 amostras
z = unifor_GenRan(10000,-10,10); %A variável z recebe valores gerados pela função uniforme 2
x = sort(z); %Atribui à variável x o vetor ordenado de z

for i = 1:k; %Laço que percorre do primeiro vetor ao vetor 10000
   Y(i) = alfa*tan(x(i)*pi)-(pi/2); %Função inversa da Cauchy 
%Utilizada para gerar as sequências com a distribuição
end
figure 1;
hist (Y,50000); %Histograma de distribuição cauchy por 50 000 divisões
title ('DISTRIBUICAO DE CAUCHY')
med = 0; %Atribui valor inicial 0 a média

for j = 1:k %Percorre o vetor de 1 a k
    med = med + Y(j); %Realiza o somatório de todos os valores encontrados
end

med = med/10000 %Divide pelo valor total
varig = 0; %Atribui valor inicial 0 a variância

for l = 1:k %Percorre o vetor de 1 a k
    varig = varig + (Y(l)-med)^2; %Calcula a variância parcial
end
varig = varig/10000 %Divide pelo número de termos para obter a variância final
axis([-200 200]) %Restringe o gráfico aos intervalos -200 e 200

for n = 1:k 
FDACacau(n) = (1/pi)*atan(x(n)/alfa) + 0.5; %FDA de Cauchy, obtida pela integração de sua FDP
end
figure 2;
plot(x, FDACacau); %Plotagem de gráfico de FDA por x
