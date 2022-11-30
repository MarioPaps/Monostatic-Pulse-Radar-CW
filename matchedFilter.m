%%Marios Papadopoulos, 4th year EEE, 01527402
%%function to implement matched filter for received signal
%parameters
%rxsig: received signal vector at point Z
% pulse: pulse compression vector

function [mf_out] = matchedFilter(rxsig,pulse)

    h_t=flip(pulse);
    mf_out= (1/length(h_t))*conv(rxsig,h_t);
    mf_out= mf_out(1:end-6); %truncate to maintain signal length as 1x11200
   
end