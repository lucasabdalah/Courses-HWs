function symb= mapping_MQAM(M,bits)


    switch M
    
    % Se estivermos trabalhando com 4-QAM
    case 4 
        switch bits
            case [0000]
                symb = 
             
        end
    
    % Se estivermos trabalhando com 16-QAM
    case 16
    
    % Se estivermos trabalhando com 64-QAM
    case 64
    
    end


% To get Gray-code back to binary it needs only to sort()