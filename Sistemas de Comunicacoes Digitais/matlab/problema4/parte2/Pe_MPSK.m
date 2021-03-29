function erro = Pe_MPSK(M,Es,N0)
    % erro = Pe_MPSK(M,Es,N0)
    % Computa a probabilidade do erro da constelacao MQAM
    %
    % SYNTAX: erro = Pe_MPSK(M,SNR_linear);
    %
    % INPUTS: 
    %       M : Quantidade de simbolos da constelacao 
    %       SNR_linear : Relacao sinal-ruido
    % OUTPUTS:
    %       erro : Erro
    %
    % Referencia: Pagina 170 - Cecilio (1a ed.)
    % Eq:    P(e) = 4 \left(1-\frac{1}{\sqrt{M}}\right) Q\left(\sqrt{\frac{3}{M-1}\frac{E_s}{N_0}}\right) - 4\left(1-\frac{1}{\sqrt{M}}\right)^2 Q^2\left(\sqrt{\frac{3}{M-1}\frac{E_s}{N_0}}\right)
    %
    %HISTORY:
    % 2021/03/27: - Lucas Abdalah.
    %
    
    erro = 2*qfunc(sqrt((2*Es)/N0) *sin(pi/M));

    end