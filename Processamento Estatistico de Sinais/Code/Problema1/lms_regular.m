function [y,e,w] = lms_regular(input)
% lms_regular - Compute the classical LMS algorithm
%
% Syntax: [y,e,w] = lms_regular(input)
%
% Long description
    
end

%% 
equalizer_length = 35;
mu = 0.4;

% Equalizer initialization
w  = zeros(equalizer_length,1); % equalizer coefficients (column) 
u  = zeros(1,equalizer_length); 
e = zeros(1,N); % error vector
num_errors=0;

%  Adaptive Equalization

% Training mode
for i = 1:training+Delta
   u  = [r(i) u(1:equalizer_length-1)];
   hat_s(i) = u*w;
   d(i) = train_data(i); % training
   e(i) = d(i)- hat_s(i);
   w = w + mu*e(i)*u';
end


def lmsPred(x,l,u,N):

    xd= np.block([np.zeros((1,l)), x]).T
    y=np.zeros((len(xd),1))
    xn=np.zeros((N+1,1))
    xn = np.matrix(xn)
    wn=np.random.rand(N+1,1)/10
    M=len(xd)

    for n in range(0,M):
        xn = np.block([[xd[n]], [xn[0:N]]]);
        y[n]= np.matmul(wn.T, xn);
        if(n>M-l-1):
            e =0;
        else:
            e=int(x[n]-y[n]);
        wn = wn + 2*u*e*xn;
        
    return y,wn;