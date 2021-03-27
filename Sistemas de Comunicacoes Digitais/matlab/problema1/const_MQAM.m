function const_QAM = const_MQAM(M,d)
    % Computa os simbolos da constelacao M-QAM, sabendo a quantidade de simbolos (M) e a distancia entre os simbolos (d).
    %
    % SYNTAX: const_QAM = const_MQAM(M,d);
    %
    % INPUTS: 
    %       M : Quantidade de simbolos da constelacao 
    %       d : Distancia entre os simbolos
    %
    % OUTPUTS:
    %       const_QAM : Simbolos da constelacao M-QAM ordenados para codigo gray
    %
    % Referencia: Pagina 58 - Cecilio (1a ed.) 
    % Componentes A_i e B_i assumem valores discretos: $\{ (2i -\sqrt{M} - 1)d \}$ 
    % 
    %HISTORY:
    % 2021/03/25: - Lucas Abdalah.
    %
    
    ai = zeros(1,sqrt(M));
    for ii=1:sqrt(M)
        ai(1,ii) = (2*ii - sqrt(M) - 1)*d;
    end
    
    const_QAM = zeros(1,M); % - Vetor que recebe o alfatbeto QAM.
    k = 1;
    
    for ii = 1:sqrt(M)      % bi - moving on the imaginary axis  
        for jj = 1:sqrt(M)  % ai - moving son the real axis  
            const_QAM(1,k) = ai(1,jj) + 1j*ai(1,ii);
            k=k+1;
        end
    end

    % This part is to display the symbol's bits
    % To make this part work is necessary to organize the M-QAM
    % modulation to sort the binary constellation and convert it into encoding
    mat_QAM = reshape(const_QAM,[sqrt(M),sqrt(M)])';
    mat_QAM(2:2:end,:) = fliplr(mat_QAM(2:2:end,:));
    aux = mat_QAM.';
    const_QAM = flip(aux(:)');

end % function