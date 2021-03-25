function [symb,bits]= demapping_MQAM(r,d)

    % M-QAM demapping
    if r <=0
        symb= -1;
        bits= 0;
    end
    if r>0
        symb = 1;
        bits= 1;
    end
    