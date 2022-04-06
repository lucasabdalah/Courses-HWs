%   Função que transforma as coordenadas esféricas (r,f,t) para coordenadas retangulares (x,y,z)
%   Função desenvolvida pelo prof. Sérgio Antenor de Carvalho 
%
%   Inicio  24/02/2010 - versão Scilab
%   Revisão 28/08/2017 - versão Matlab
%
%   Dados de entrada: coordenadas esféricas
%   r -> coordenada r; f -> coordenada fi em radianos está entre 0 e 2pi; t -> coordenada teta em radianos está entre 0 e pi
%
%   Dados de saída: coordenadas retangulares
%   x -> coordenada x; y -> coordenada y; z -> coordenada z

function [x,y,z] = CooEsfRet(r,t,f)
   x = r .* sin(t) .* cos(f); 
   y = r .* sin(t) .* sin(f);
   z = r .* cos(t);
end % função CooEsfRet