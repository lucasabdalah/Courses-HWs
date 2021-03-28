function erro = Pe_MQAM(M,Es,N0)
% erro = Pe_MQAM(M,Es,N0)
% Computa a probabilidade do erro da constelacao MQAM
%
% SYNTAX: erro = Pe_MQAM(M,SNR_linear);
%
% INPUTS: 
%       M : Quantidade de simbolos da constelacao 
%       SNR_linear : Relacao sinal-ruido
% OUTPUTS:
%       erro : Erro
%
% Referencia: Pagina 143 e 144 - Cecilio (1a ed.)
% Eq:    P(e) = 4 \left(1-\frac{1}{\sqrt{M}}\right) Q\left(\sqrt{\frac{3}{M-1}\frac{E_s}{N_0}}\right) - 4\left(1-\frac{1}{\sqrt{M}}\right)^2 Q^2\left(\sqrt{\frac{3}{M-1}\frac{E_s}{N_0}}\right)
%
%HISTORY:
% 2021/03/26: - Lucas Abdalah.
%

erro = 4*(1-(1/sqrt(M)))*first_part(M,Es,N0)-(4*(1-(1/sqrt(M)))^2)*second_part(M,Es,N0);
end

%% To make simpler the expression
function primeiro_termo = first_part(M,Es,N0)
    primeiro_termo = qfunc(sqrt( (3*Es)/(N0*(M-1)) ));
end
function segundo_termo = second_part(M,Es,N0)
    segundo_termo = (first_part(M,Es,N0))^2;
    %% Teste da equcao simplificada, com SNR alta
    % segundo_termo = 0;
end