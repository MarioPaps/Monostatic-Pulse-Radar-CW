%% load data
load("BackscatterData.mat");
backscatmat=cell2mat(BackscatterData);
%% setup
A=150;
Tc=28e-9;
Fc=15e9; %Ku band ranges from 12 to 18 GHz, so the middle is 15 GHz
c=3e8;
lambda=c/Fc;
d=lambda/2;
rx=(-22:1:22)*d;
pulsecomp=[-1,-1,-1,+1,+1,-1,+1];
Fn=3.1068;
kB=1.28e-23;
To=290; %in Kelvin
B=1/Tc; %pulse compression bandwidth
Pn= kB*To*Fn*B;
%store the backscatter data for each angle
const=length(BackscatterData);
backscatter_angle= zeros(45,11200,length(BackscatterData));
for ind=1:const
    backscatter_angle(:,:,ind)= cell2mat(BackscatterData(ind));
end
%% transmit noise in another direction to determine threshold
deg=30; 
thetasteering= deg2rad(deg); 
Txout= zeros(length(rx),const); %do not transmit anything since the branches will be skipped
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
threshold=specifyThreshold(Rxout,7);
%% angle estimation
%sweep through the angles to identify maxima
max_storer=zeros(const,1);
for theta=30:150
   thetast=deg2rad(theta);
   [Rxout,~]=Rxfun(backscatter_angle(:,:,theta-29),thetast,rx,lambda);

   %matched filter
   mfout=matchedFilter(Rxout,pulsecomp);
   %apply threshold
   mfout(abs(mfout)<threshold)=0;

   %average across 8 PRIs  
   pri_sigs= reshape(abs(mfout),[],8);
   val_mean= mean(pri_sigs,2);
   max_storer(theta-29)= max(val_mean);
end
figure;
plot((30:150),max_storer);
xlabel('Angle(degrees)'); ylabel('Magnitude(V)'); title('Angle Sweep Plot');
%% identify number of targets
%a maximum is found at angle=137 degrees so the backscatter data will be plotted for this angle
[~,deg]= max(max_storer);
deg=deg+29;
thetast= deg2rad(deg);
[Rxout,~]=Rxfun(backscatter_angle(:,:,deg-29),thetast,rx,lambda);

%plot signal at point Z
figure;
semilogy(abs(Rxout));
xlabel('Time(Tc)'); ylabel('Magnitude(V)'); title('Received Signal at Point Z');
points=(1:8)*1400;
xline(points,'--r');
yline(threshold,'-','Threshold');

%matched filter
mfout=matchedFilter(Rxout,pulsecomp);
%apply threshold
mfout(abs(mfout)<threshold)=0;

%average across 8 PRIs  
pri_sigs= reshape(abs(mfout),[],8);
val_mean= mean(pri_sigs,2);
figure;
plot((1:length(Rxout)/8),val_mean);
xlabel('Time(Tc)'); ylabel('Magnitude(V)'); title('Signal Average Across PRIs');
%% Parameter estimation for the 4 targets
num=4;
[peak,delays]= maxk(val_mean,num);
techo= (delays-7).*Tc;
R_estimated= (c.*techo)./2;
GTx=45;
GRx=45;
sqrt_RCS_est= (peak.*sqrt((4*pi)^3).* R_estimated.^2) ./ (A*GTx*GRx*lambda);
RCS_estimated= sqrt_RCS_est.^2;




