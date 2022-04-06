%   Função que transforma um vetor em coordenadas cilíndricas (ro,fc,zc) para um vetor em coordenadas retangulares (x,y,z)
%   Função desenvolvida pelo prof. Sérgio Antenor de Carvalho 
%
%   Inicio  24/02/2010 - versão Scilab
%   Revisão 28/08/2017 - versão Matlab
%
%   Dados de entrada: vetor em coordenadas cilíndricas
%   ro, fc, z     -> coordenadas ro, f, z do ponto onde o vetor é definido em coordenadas cilíndricas
%   Aro, Afc, Azc -> componentes do vetor no sistema de coordenadas esféricas 
%
%   Dados de saída: vetor em coordenadas retangulares
%   x, y, z    -> coordenadas x, y, z do ponto onde o vetor é definido
%   Ax, Ay, Az -> componentes do vetor no sistema de coordenadas retangulares  

function [x,y,z,Ax,Ay,Az] = VetCilRet(ro,fc,zc,Aro,Afc,Azc)
   [x,y,z] = CooCilRet(ro,fc,zc);  % calcula as coordenadas retangulares do ponto onde o vetor é definido 
   senf = sin(fc); cosf = cos(fc); % cálculos auxiliares
   Ax = Aro.* cosf - Afc * senf;
   Ay = Aro * senf + Afc * cosf;
   Az = Azc;
end % função VetCilRet
