function [symb,bits]= pam4_demapping(r,d)

% 4-PAM demapping
if r <=0
    symb= -1;
    bits= 0;
end
if r>0
    symb = 1;
    bits= 1;
end
