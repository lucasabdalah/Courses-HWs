%   Função que transforma as coordenadas cilíndricas (r,f,z) para coordenadas retangulares (x,y,z) 
%   Função desenvolvida pelo prof. Sérgio Antenor de Carvalho 
%
%   Inicio  24/02/2010 - versão Scilab
%   Revisão 28/08/2017 - versão Matlab
%
%   Dados de entrada: coordenadas cilíndricas
%   ro -> coordenada ro; fc -> coordenada fi em radianos está entre 0 e 2pi; zc -> coordenada z
%
%   Dados de saída: coordenadas retangulares
%   x -> coordenada x; y -> coordenada y; z -> coordenada z

function [x,y,z] = CooCilRet(ro,fc,zc)
   x = ro.* cos(fc); 
   y = ro.* sin(fc);
   z = zc;
end % função CooCilRet