%%Marios Papadopoulos, 4th year EEE, 01527402
%%function to generate backscatter data for given target direction
%parameters
% in: input signal
% targ_theta, targ_range, targ_RCS, targ_type : target direction, range,
% average RCS and type vectors
%rx: x-coordinate vector of antenna array
% lambda,c,Fc: wave length, velocity of light, and carrier frequency
% Tc: chip period
% Pn: noise power
function[backs_out,Tx_man,Rx_man]=backscatterfn(in,targ_theta,targ_range,targ_RCS,targ_type,rx,lambda,c,Fc,Tc,Pn)
    GTx=1;
    GRx=1;
    L=length(targ_theta);
    acc=zeros(size(in));
    if(L~=0)
        for ind=1:L
            psicurr= computepsi(rx,lambda,targ_theta(ind));
            Tx_man= exp(1i*psicurr);
            Rx_man= exp(-1i*psicurr);
            
            sig_Tx= transpose(Tx_man)*in;
            techo= 2*targ_range(ind)/c; %delay in chips is translated to ceil(techo/Tc)
            delay= ceil(techo/Tc);
            sig_del=[zeros(1,delay),sig_Tx(1:end-delay)];

            switch targ_type(ind)
                case "simple"
                    magn_beta=sqrt(ones(1,length(in)));
                case "similar"
                    temp_beta= fSwerling12rnd(targ_RCS(ind),8,'Amplitude');
                    magn_beta= repelem(temp_beta,length(in)/length(temp_beta));
                case "large"
                    temp_beta= fSwerling34rnd(targ_RCS(ind),8,'Amplitude');
                    magn_beta= repelem(temp_beta,length(in)/length(temp_beta));
            end
            random_phase= 2*pi* rand(1,8);
            random_phase= repelem(random_phase,length(in)/8);

beta= sqrt((GTx*GRx)/(4*pi)^3) *(lambda/(targ_range(ind)^2)).*magn_beta.*exp(-1i*2*pi*Fc*(2*targ_range(ind)/c)) .*exp(1i*random_phase);
            sig_att= beta.*sig_del;
            sig_rx=  Rx_man * sig_att;
            acc= acc + sig_rx;
        end
    end

        %define noise 
         noise= sqrt(Pn/2)*( randn(size(in)) + 1i*randn(size(in)) );
        
         %the first 7 Tc of every PRI must be 0
         pos_zeros=[];
         for ind=1:8
            newvec= 1+(ind-1)*1400:1: 7+(ind-1)*1400;
            pos_zeros=[pos_zeros, newvec];
         end

         backs_out= acc+noise; 
         backs_out(1:45,pos_zeros)=0; %assign zeros
end









