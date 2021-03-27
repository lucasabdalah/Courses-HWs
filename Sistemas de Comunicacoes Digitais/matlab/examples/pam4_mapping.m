function symb= pam4_mapping(bits,d)

% 4-PAM mapping
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