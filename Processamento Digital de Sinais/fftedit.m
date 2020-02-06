function [filefft, freq] =  fftedit(file, fs)
    % Run the Fast Fourier Transform and do the abs of the complex numbers
  
    % To make it faster and time axes 
    n = length(file); 
    t = n/fs; %File length and time range
    k = 0:n-1; 
    freq=k/t; % Creates the frequency range
    
    % Help to cut filefft and frequency vectors
    periodo=ceil(n/2); 
    freq=freq(1:periodo); 
  
    %Do the abs of fftn data and cut
    filefft = abs(fft(file));
    filefft=filefft(1:periodo);
end