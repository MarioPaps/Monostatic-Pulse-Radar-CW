%%Marios Papadopoulos, 4th year EEE, 01527402
%%function to compute phase shifter vector for given target direction
%parameters
%rx: x-coordinate vector of antenna array
% lambde: wavelength
% thetast: steering direction
function[psi]= computepsi(rx,lambda, thetast)
   k=((2*pi)/lambda) * cos(thetast);
   psi= rx'*k;
 
end