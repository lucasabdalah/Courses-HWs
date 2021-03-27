% Problema 1 - Parte 4
function [symb,bits]=demapping_MQAM(r,M,d)
% Computa o simbolo e os bits da mensagem de acordo com a r (coordenada recebida), M (Quantidade de simbolos da constelacao) e distancia da constelacao (d).
%
% SYNTAX: [symb,bits]= demapping_MQAM(r,d,M);
%
% INPUTS: 
%       r : Coordenada recebida
%       M : Quantidade de simbolos da constelacao 
%       d : distancia da constelacao
% OUTPUTS:
%       symb : Simbolo da constelacao
%       bits : Simbolos binarios com codificao de Gray
%
% EXAMPLE:
%       r = 0.5+1j*0.5;
%       M = 4;
%       d = sqrt(2)/2;
%       [symb,bits]=demapping_MQAM(r,M,d);
%
%HISTORY:
% 2021/03/26: - Lucas Abdalah.
%

const_QAM = const_MQAM(M,d);    % Constellation M-QAM
d_E = zeros(1,M);
for ii = 1:M
    r0 = const_QAM(ii);
    d_E(ii) = sqrt((real(r)-real(r0))^2 + (imag(r)-imag(r0))^2);
end

%% Get the minimum Euclidean distance for all symbols
[r0,get_position] = min(d_E);
gray_alfabeto = gray_const(M,false);    % Gray alphabet

%% Returning the values
symb = const_QAM(get_position);
bits = gray_alfabeto(get_position,:);

end