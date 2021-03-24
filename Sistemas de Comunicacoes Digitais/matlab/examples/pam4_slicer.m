function [symb_dec,bits_dec]= pam4_slicer(bits,d)
    
    % Signal lenght
    L=length(bits)
    bits_dec = zeros(1,L);
    % Dividir o sinal em 4
    for windowing = 1:4:L
        
        bits(windowing);


    
    % 4-PAM slicing
    if bits==[0 0]
        symb= -3*d;
    end
    if bits==[0 1]
        symb= -1*d;
    end
    if bits==[1 1]
        symb = 1*d;
    end
    if bits==[1 0]
        symb= 3*d;
    end