function d = d_MQAM(M,E_media)
    % Computa a distancia entre os simbolos da constelacao M-QAM, sabendo a quantidade de simbolo M e a energia media (E_media).
    %
    % SYNTAX: d = d_MQAM(M,E_media);
    %
    % INPUTS: 
    %       M : Quantidade de simbolos da constelacao 
    %       E_media : Energia media da constelacao
    %
    % OUTPUTS:
    %       d : Distancia entre os simbolos
    %
    % Referencia: Pagina 142, eq.3.186 - Cecilio (1a ed.)
    % $\sqrt{\frac{3 \mathcal{E}_{media}}{2(M-1)}} $
    %
    % HISTORY:
    % 2021/03/24: - Lucas Abdalah.
    %
    
    d = sqrt(3*E_media/(2*(M-1)));
    
end