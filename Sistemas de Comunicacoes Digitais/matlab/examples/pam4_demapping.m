function [symb,bits]= pam4_demapping(r,d)

% 4-PAM demapping
if r <= -2*d
    symb= -3*d;
    bits= [0 0];
end
if r>-2*d && r<= 0
    symb = -1*d;
    bits= [0 1];
end
if r>0 && r<= 2*d
    symb = 1*d;
    bits= [1 1];
end
if r > 2*d
    symb= 3*d;
    bits= [1 0];
end