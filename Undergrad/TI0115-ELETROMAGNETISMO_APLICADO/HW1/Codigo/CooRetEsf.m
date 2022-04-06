%   Função que transforma as coordenadas retangulares (x,y,z) para coordenadas esféricas (r,f,t)
%   dados de entrada: coordenadas retangulares
%   Função desenvolvida pelo prof. Sérgio Antenor de Carvalho 
%
%   Inicio  24/02/2010 - versão Scilab
%   Revisão 28/08/2017 - versão Matlab
%
%   Dados de entrada
%   x -> coordenada x; y -> coordenada y; z -> coordenada z
%
%   Dados de saída: coordenadas esféricas
%   r -> coordenada r; f -> coordenada fi em radianos está entre 0 e 2pi; t -> coordenada teta em radianos está entre 0 e pi

function [r,t,f] = CooRetEsf(x,y,z)
   r = sqrt(x.^2 + y.^2 + z.^2); 
   t = acos(z./r);
   f = FI(x,y);
end% função CooRetEsf