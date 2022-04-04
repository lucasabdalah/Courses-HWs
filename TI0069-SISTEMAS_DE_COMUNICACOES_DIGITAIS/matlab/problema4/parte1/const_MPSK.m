function const_PSK = const_MPSK(M,E_s)
% Computa os simbolos da constelacao M-PSK, sabendo a quantidade de simbolos (M) e a distancia entre os simbolos (d).
%
% SYNTAX: const_PSK = const_MPSK(M,d);
%
% INPUTS: 
%       M : Quantidade de simbolos da constelacao 
%       d : Distancia entre os simbolos
%
% OUTPUTS:
%       const_PSK : Simbolos da constelacao M-QAM ordenados para codigo gray
%
% Referencia: Pagina 150 - Cecilio (1a ed.) 
%
% 
%HISTORY:
% 2021/03/27: - Lucas Abdalah.
%

for ii=1:(M/2)
    phi_i(1,ii) = ((2*ii - 1)*pi)/M;
end

const_PSK = zeros(1,M); % - Vetor que recebe o alfatbeto QAM.

k = 1;
for ii = 1:(M/2)      % bi - moving on the imaginary axis  
    const_PSK(1,ii) = (+cos(phi_i(1,ii)) + 1j*sin(phi_i(1,ii)));
end
const_PSK(1,(M/2+1):end) = -const_PSK(1,1:(M/2));
const_PSK = sqrt(E_s)*(const_PSK);

end % function