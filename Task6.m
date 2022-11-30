%Fourth scan-targets 1,2,3  present
addpath("Yunhao Functions");
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
%% signal transmission
Rx_total=[];
max_storer=zeros(1,121);
for deg=30:150 %30:150
    thetasteering=deg2rad(deg); 
    [Txout,Txmanifold]= Txfun(in_dwell,thetasteering,rx,lambda);
   
   %target parameters
    targ_thetas=[40,70,120]; %in deg
    targ_thetas=deg2rad(targ_thetas);
    targ_range=[2e3,3e3,2.5e3];
    targ_RCS=[1,5,4.5];
    targ_type=["simple","similar","large"];

    %the signal is then sent to the backscatter
    [backsc_out,target_tx,target_rx]= backscatterfn(Txout,targ_thetas,targ_range,targ_RCS,targ_type,rx,lambda,c,Fc,Tc,Pn);

    %the backscatter output is sent to the Rx
    [Rxout,Rxmanifold]= Rxfun(backsc_out,thetasteering,rx,lambda);

    %average across 8 PRIs
    pri_sigs= reshape(abs(Rxout),[],8);
    val_mean= mean(pri_sigs,2);
    max_storer(deg-29)=max(val_mean);
end
%plot angle sweep
plot((30:150),max_storer);
xlabel('Angle(degrees)'); ylabel('Magnitude(V)'); title('Angle Sweep Plot');

[~,angle_est]= maxk(max_storer,3); %find k maxima
angle_est= angle_est+29;

%vectors to store estimated target parameters
techo=zeros(3,1);
R_estimated=zeros(3,1);
peak= zeros(3,1);
RCS_estimated= zeros(3,1);
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
threshold=specifyThreshold(Rxout,6);
%% Range and RCS estimation for 40 degrees
deg=40;
thetasteering=deg2rad(deg); 
[Txout,Txmanifold]= Txfun(in_dwell,thetasteering,rx,lambda);
targ_thetas=[40];
targ_thetas=deg2rad(targ_thetas);
targ_range=2e3;
targ_RCS=1;
targ_type="simple";

%the signal is then sent to the backscatter
[backsc_out,target_tx,target_rx]= backscatterfn(Txout,targ_thetas,targ_range,targ_RCS,targ_type,rx,lambda,c,Fc,Tc,Pn);

%the backscatter output is sent to the Rx
[Rxout,Rxmanifold]= Rxfun(backsc_out,thetasteering,rx,lambda);

%plot at point Z
figure;
semilogy(abs(Rxout));
xlabel('Time(Tc)'); ylabel('Magnitude(V)'); title('Received Signal at Point Z');
points=(1:8)*length(in_PRI);
xline(points,'--r');

%matched filter and threshold
mfout= matchedFilter(Rxout,pulsecomp);
mfout(abs(mfout)<threshold)=0;

%average across 8 PRIs  
pri_sigs= reshape(abs(mfout),[],8);
val_mean= mean(pri_sigs,2);
figure;
plot((1:1400),val_mean);
xlabel('Time Index(Tc)'); ylabel('Magnitude(V)'); title('Signal Average Across PRIs');
[peak(1),pos_delay]= max(val_mean);
techo(1)= (pos_delay-7)*Tc; %the first 7 samples from the convolution are 0
R_estimated(1)= c*techo(1)/2;  
%% Range and RCS estimation for 70 degrees
deg=70;
thetasteering=deg2rad(deg); 
[Txout,Txmanifold]= Txfun(in_dwell,thetasteering,rx,lambda);
targ_thetas=70;
targ_thetas=deg2rad(targ_thetas);
targ_range=3e3;
targ_RCS=5;
targ_type="similar";

%the signal is then sent to the backscatter
[backsc_out,target_tx,target_rx]= backscatterfn(Txout,targ_thetas,targ_range,targ_RCS,targ_type,rx,lambda,c,Fc,Tc,Pn);

%the backscatter output is sent to the Rx
[Rxout,Rxmanifold]= Rxfun(backsc_out,thetasteering,rx,lambda);

%plot at point Z
figure;
semilogy(abs(Rxout));
xlabel('Time(Tc)'); ylabel('Magnitude(V)'); title('Received Signal at Point Z');
points=(1:8)*length(in_PRI);
xline(points,'--r');

%matched filter and threshold
mfout= matchedFilter(Rxout,pulsecomp);
mfout(abs(mfout)<threshold)=0;

%average across 8 PRIs  
pri_sigs= reshape(abs(mfout),[],8);
val_mean= mean(pri_sigs,2);
figure;
plot((1:1400),val_mean);
xlabel('Time Index(Tc)'); ylabel('Magnitude(V)'); title('Signal Average Across PRIs');
[peak(2),pos_delay]= max(val_mean);
techo(2)= (pos_delay-7)*Tc; %the first 7 samples from the convolution are 0
R_estimated(2)= c*techo(2)/2;
%% Range and RCS estimation for 120 degrees
deg=120;
thetasteering=deg2rad(deg); 
[Txout,Txmanifold]= Txfun(in_dwell,thetasteering,rx,lambda);
targ_thetas=120;
targ_thetas=deg2rad(targ_thetas);
targ_range=2.5e3;
targ_RCS=4.5;
targ_type="large";

%the signal is then sent to the backscatter
[backsc_out,target_tx,target_rx]= backscatterfn(Txout,targ_thetas,targ_range,targ_RCS,targ_type,rx,lambda,c,Fc,Tc,Pn);

%the backscatter output is sent to the Rx
[Rxout,Rxmanifold]= Rxfun(backsc_out,thetasteering,rx,lambda);

%plot at point Z
figure;
semilogy(abs(Rxout));
xlabel('Time(Tc)'); ylabel('Magnitude(V)'); title('Received Signal at Point Z');
points=(1:8)*length(in_PRI);
xline(points,'--r');

%matched filter and threshold
mfout= matchedFilter(Rxout,pulsecomp);
mfout(abs(mfout)<threshold)=0;

%average across 8 PRIs  
pri_sigs= reshape(abs(mfout),[],8);
val_mean= mean(pri_sigs,2);
figure;
plot((1:1400),val_mean);
xlabel('Time Index(Tc)'); ylabel('Magnitude(V)'); title('Signal Average Across PRIs');
[peak(3),pos_delay]= max(val_mean);
techo(3)= (pos_delay-7)*Tc; %the first 7 samples from the convolution are 0
R_estimated(3)= c*techo(3)/2;
%% RCS estimation
GTx=45;
GRx=45;
sqrt_RCS_est= (peak.*sqrt((4*pi)^3).* R_estimated.^2) ./ (A*GTx*GRx*lambda);
RCS_estimated= sqrt_RCS_est.^2;
%% Display results
disp('Estimated t_echo are:');
disp(techo);
disp('Estimated ranges are:');
disp(R_estimated);
disp('Estimated RCS are:');
disp(RCS_estimated);

