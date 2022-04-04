%   Função que transforma um vetor em coordenadas esféricas (r,f,t) para um vetor em coordenadas retangulares (x,y,z)
%   Função desenvolvida pelo prof. Sérgio Antenor de Carvalho 
%
%   Inicio  24/02/2010 - versão Scilab
%   Revisão 28/08/2017 - versão Matlab
%
%   Dados de entrada: vetor em coordenadas esféricas
%   r, t, f    -> coordenadas r, t, f do ponto onde o vetor é definido em coordenadas esféricas
%   Ar, At, Af -> componentes do vetor no sistema de coordenadas esféricas 
%
%   Dados de saída vetor em coordenadas retangulares
%   x, y, z    -> coordenadas x, y, z do ponto onde o vetor é definido
%   Ax, Ay, Az -> componentes do vetor no sistema de coordenadas retangulares  

function [x,y,z,Ax,Ay,Az] = VetEsfRet(r,t,f,Ar,At,Af)
   [x,y,z] = CooEsfRet(r,t,f); % calcula as coordenadas retangulares do ponto onde o vetor é definido 
   sent = sin(t); cost = cos(t); senf = sin(f); cosf = cos(f); % cálculos auxiliares
   Ax = Ar .* sent .* cosf + At .* cost .* cosf - Af .* senf;       % componente Ax
   Ay = Ar .* sent .* senf + At .* cost .* senf + Af .* cosf;       % componente Ay
   Az = Ar .* cost - At .* sent;                                 % componente Az
end % função VetEsfRet