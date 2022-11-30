%%Marios Papadopoulos, 4th year EEE, 01527402
%%function to implement receiver phase shifter and weights
%parameters
%insig: received signal vector after baseband processing at Rx
%thetast: steering angle in radians
%rx: x-coordinate vector of antenna array
%lambda: wavelength

function[Rx_out,Rxmanifold]= Rxfun(insig, thetast,rx,lambda)
    
    psicurr= computepsi(rx,lambda,thetast);
    Rxmanifold= exp(-1i.*psicurr);
    Rxmanifold= conj(Rxmanifold);
    
    inter_sig= Rxmanifold.*insig; %signal at point D
    w=ones(length(rx),1);
    
    Rx_out= w'*inter_sig;
end