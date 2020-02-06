function y = unifor_GenRan(N = 1, inf = 0, sup = 1); %Exerce a mesma função de unifor_Gen
%Alterando apenas os limites inferior e superior para 0 e 1
clc;
%Diferente do problema anterior, a exponencial acontece nos intervalos de 0 a 1
a = 7;
b = 11;
m = 2^32;
seed = clock;
seed = ((10^3)*seed(6));
x = zeros(1,N);
x(1) = seed;

for h = 2:N
   x(h) = mod((a*x(h-1)+b),m);
end 
x = x/(m-1);

y = (x*(sup-inf))+inf;

endfunction
%Utilizada para alimentar as outras funções de distribuição (Gaussiana, exponencial e Cauchy)