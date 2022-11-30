%First scan-only noise is sent
addpath("Yunhao Functions\");
A=150;
Tp=196e-9;
Tc=28e-9;
B=1/Tc; %pulse compression bandwidth
pulsecomp=[-1,-1,-1,+1,+1,-1,+1];
M=199;
Pfa=1e-3;

Fc=15e9; %Ku band ranges from 12 to 18 GHz, so the middle is 15 GHz
c=3e8;
lambda=c/Fc;
Fn=3.1068;
kB=1.28e-23;
To=290; %in Kelvin
Pn= kB*To*Fn*B; %noise power

d=lambda/2;
rx=(-22:1:22)*d;
ry=zeros(1,length(rx));
rz=zeros(1,length(rx));
r=[rx;ry;rz]; %same antenna array at Tx and Rx
in= A.*pulsecomp; %it contains 7Tc ->1 cell=1 Tc
in_PRI=[in,zeros(1,M*length(pulsecomp))]; % 1PRI 1x1400
in_dwell= repmat(in_PRI,[1 8]); %8 PRIs => 1x11200
%% 1st scan- no targets

%signal is steered towards a random direction
deg=30; 
thetasteering= deg2rad(deg); 
Txout= Txfun(in_dwell,thetasteering,rx,lambda);

%target parameters
targ_thetas=[];
targ_range=[];
targ_RCS=[];
targ_type=[];

%the signal is then sent to the backscatter
backsc_out= backscatterfn(Txout,targ_thetas,targ_range,targ_RCS,targ_type,rx,lambda,c,Fc,Tc,Pn);

%the backscatter output is sent to the Rx
Rxout= Rxfun(backsc_out,thetasteering,rx,lambda);

%plot noise samples magnitude at point Z
figure;
semilogy(abs(Rxout));
xlabel('Time Index (Tc)'); ylabel('Magnitude(V)'); title('Noise Signal at Point Z');
points=(1:8)*length(in_PRI);
xline(points,'--r');

%estimate and plot pdf of received samples
[N,edges]= histcounts(abs(Rxout),'Normalization','pdf');
figure;
histogram(abs(Rxout),length(N),'Normalization','pdf');
xlabel('Magnitude(V)'); ylabel('Probability'); title('Noise PDF Estimate');
hold on;
%estimate Pn at point Z
Pn_est= (1/length(Rxout)) * (Rxout)*ctranspose(Rxout);
Pn_est= real(Pn_est);

%noise magnitude samples have a Rayleigh pdf
p_theor= raylpdf(abs(Rxout),sqrt(Pn_est/2));
plot(abs(Rxout),p_theor,'o','LineStyle', 'none','Color','g');
legend('Empirical PDF','Theoretical Rayleigh PDF');
hold off;

%% 
%compute threshold
threshold=specifyThreshold(Rxout,3)

