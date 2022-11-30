%%Marios Papadopoulos, 4th year EEE, 01527402
%%function to implement receiver phase shifter and weights
%parameters
%insig: received signal vector after baseband processing at Rx
%thetast: steering angle in radians
%rx: x-coordinate vector of antenna array
%lambda: wavelength

function[Tx_out,Txmanifold]= Txfun(insig,thetast,rx,lambda)
    w=ones(length(rx),1);
    sig_inter= w*insig; %apply weights ->45x11200
    psicurr= computepsi(rx,lambda,thetast);
    
    Txmanifold= exp(1i.*psicurr);
    Txmanifold= conj(Txmanifold);
    Tx_out= Txmanifold.* sig_inter; %apply Tx manifold - 45x11200 signal out
    
end