%define constants
Tc=28e-9;
Fc=15e9; %Ku band ranges from 12 to 18 GHz, so the middle is 15 GHz
c=3e8;
Tp=196e-9;
M=199; %range bins
PRI= (M+1)* Tp;
Tdwell= 8*PRI;

%define antenna array
lambda=c/Fc;
d=lambda/2;
rx=(-22:1:22)*d;
ry=zeros(1,length(rx));
rz=zeros(1,length(rx));
r=[rx;ry;rz]; %same antenna array at Tx and Rx

%compute steering vectors
thetamain=[40,70,120]; %in degrees
thetamain= deg2rad(thetamain);
phimain=0;
psimain=zeros(length(rx),3);

%psimain is the vector of phase shifters
for iter=1:length(thetamain)
    psimain(:,iter)= (2*pi/lambda)*rx * cos(thetamain(iter));
end

psimain= mod(psimain,2*pi); %vector of phase shifters expressed in the [0,2*pi] range
%psimain= rad2deg(psimain); %vector of phase shifters expressed in the [0,360] range in degrees
%code to print phase shifter matrix to a file
% file_ID = fopen('testfile.txt','w');
% fprintf(file_ID,'%.4f \n %.4f \n %.4f \n',psimain);
% fclose(file_ID);
% type('testfile.txt')
%%  plot array patterns
for iter=1:3
    Txmanifoldmain= exp(1i.*psimain(:,iter));
    Rxmanifoldmain= exp(-1i.*psimain(:,iter));
    thetas=(0:360);
    thetas= deg2rad(thetas);
    psicurr= computepsi(rx,lambda,thetas);
    Tx_curr= exp(1i.*psicurr);
    Rx_curr= exp(-1i.*psicurr);
    gain= Rxmanifoldmain' * Rx_curr;
    gain= abs(gain);
    figure;
    plot((0:180),gain(1:181));
    xlabel('Azimuth(degrees)'); ylabel('Gain'); title('ULA Gain Pattern');
    figure;
    axis=linspace(0,2*pi,length(gain));
    polarplot(axis, gain);
end