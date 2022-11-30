%Second scan-target 1 is present
addpath("Yunhao Functions");
A=1000;
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

d=lambda/2;
rx=(-22:1:22)*d;
ry=zeros(1,length(rx));
rz=zeros(1,length(rx));
r=[rx;ry;rz]; %same antenna array at Tx and Rx

%transmitted signal
in= A.*pulsecomp; %it contains 7Tc ->1 cell=1 Tc
in_PRI=[in,zeros(1,M*length(pulsecomp))]; % 1PRI 1x1400
in_dwell= repmat(in_PRI,[1 8]); %8 PRIs => 1x11200
PTx= (1/length(in_dwell))* (in_dwell)*(in_dwell)';
%in_dwell=ones(1,11200);
Pn= kB*To*Fn*B; %noise power
%% scan over all angles in the range 30-150
Rx_total=[];
max_storer=zeros(1,121);

for deg=30:150 %30:150
    thetasteering=deg2rad(deg); 
    [Txout,Txmanifold]= Txfun(in_dwell,thetasteering,rx,lambda);
    %target parameters
    targ_thetas=40;
    targ_thetas=deg2rad(targ_thetas);
    targ_range=2e3;
    targ_RCS=1;
    targ_type="simple";
    
    %the signal is then sent to the backscatter
    [backsc_out,target_tx,target_rx]= backscatterfn(Txout,targ_thetas,targ_range,targ_RCS,targ_type,rx,lambda,c,Fc,Tc,Pn);

    %the backscatter output is sent to the Rx
    [Rxout,Rxmanifold]= Rxfun(backsc_out,thetasteering,rx,lambda);

    %check= ctranspose(Rxmanifold)* target_rx * ctranspose(target_tx)* Txmanifold
%     threshold= 4.4692e-6;
%     figure;
%     semilogy(abs(Rxout));
%     xlim([1 1400]);
%     points=(1:8)*length(in_PRI);
%     xline(points,'--r');
    %yline(threshold);
%     
%     %apply matched filter
%     mfout=matchedFilter(Rxout,pulsecomp);
%     mfout= abs(mfout);
%    % threshold= 4.4692e-6;
% %    
%     figure;
%     hold on;
%     plot(mfout);
%     points=(1:8)*length(in_PRI);
%     xline(points,'--r');
%     yline(threshold);
    
    %average across 8 PRIs
    pri_sigs= reshape(abs(Rxout),[],8);
    val_mean= mean(pri_sigs,2);
    max_storer(deg-29)=max(val_mean);
end

plot((30:150),max_storer);
xlabel('Angle(degrees)'); ylabel('Magnitude(V)'); title('Angle Sweep Plot');
% [~,angle_est]= maxk(max_storer,1); %create a function to find k maxima
% angle_est= angle_est+29;
%% transmit noise in another direction to determine threshold

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

%compute threshold
threshold=specifyThreshold(Rxout,4);

%% target angle has been estimated at 40 degrees - range and RCS estimation
angle_est=40;
thetasteering= deg2rad(angle_est); 
Txout= Txfun(in_dwell,thetasteering,rx,lambda);

%target parameters
targ_thetas=40;
targ_thetas=deg2rad(targ_thetas);
targ_range=2e3;
targ_RCS=1;
targ_type="simple";

%the signal is then sent to the backscatter
backsc_out=backscatterfn(Txout,targ_thetas,targ_range,targ_RCS,targ_type,rx,lambda,c,Fc,Tc,Pn);

%the backscatter output is sent to the Rx
[Rxout,Rxmanifold]= Rxfun(backsc_out,thetasteering,rx,lambda); %vector signal at point Z

figure;
semilogy(abs(Rxout));
xlabel('Time Index(Tc)'); ylabel('Magnitude(V)'); title('Received Signal at Point Z');
points=(1:8)*length(in_PRI);
xline(points,'--r');

PRx= (1/length(Rxout))* (Rxout)* (Rxout)';
PRx= abs(PRx);

%apply matched filter
mfout=matchedFilter(Rxout,pulsecomp);

figure;
semilogy(abs(mfout));
xlabel('Time Index(Tc)'); ylabel('Magnitude(V)'); title('Signal after Matched Filter');
points=(1:8)*length(in_PRI);
xline(points,'--r');
yline(threshold,'-','Threshold');
%apply threshold
mfout(abs(mfout)<threshold)=0;

%average across 8 PRIs  
pri_sigs= reshape(abs(mfout),[],8);
val_mean= mean(pri_sigs,2);
figure;
plot((1:1400),val_mean);
xlabel('Time Index(Tc)'); ylabel('Magnitude(V)'); title('Signal Average Across PRIs');

[peak,pos_delay]= max(val_mean);
techo= (pos_delay-7)*Tc; %the first 7 samples from the convolution are 0,leading to a shift to the right
R_estimated= c*techo/2;
%% RCS estimation
GTx=45;
GRx=45;
sqrt_RCS_est= (peak*sqrt((4*pi)^3).* R_estimated.^2) ./ (A*GTx*GRx*lambda);
RCS_estimated= sqrt_RCS_est.^2;
%% Display results
disp('Estimated t_echo are:');
disp(techo);
disp('Estimated ranges are:');
disp(R_estimated);
disp('Estimated RCS are:');
disp(RCS_estimated)