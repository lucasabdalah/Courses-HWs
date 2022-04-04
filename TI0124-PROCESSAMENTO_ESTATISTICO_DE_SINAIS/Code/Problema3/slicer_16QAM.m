function symbol = slicer_16QAM(r);
% slicer_16QAM - Compute the closest entry in a 16QAM set
%
% Syntax: y = slicer_16QAM(u);
%
% 2021/08/26 - Lucas Abdalah

%% Compute the values
if abs(real(r)) <= 2
   a = 1; 
else 
   a = 3;
end 

if abs(imag(r)) <= 2
   b = 1; 
else 
   b = 3;
end 

%% Compute the signal
if real(r) >= 0  
   a_signal = 1;   
else             
   a_signal = -1;  
end 

if imag(r) >= 0 
   b_signal = 1; 
else 
   b_signal = -1; 
end 

symbol = a*a_signal + j*b*b_signal;
