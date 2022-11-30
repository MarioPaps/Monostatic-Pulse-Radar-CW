clc;clear;close all
%% UCA Pattern
% Written by Yunhao Liu
% Last modified 08/03/2021
N = 8; % Number of array element
lambda = 0.1; % wavelength
R = 0.0653;   % UCA radius

% UCA Cartesian coordinates centered at origin
r_UCA = [R,R*cos(pi/4),0,R*cos(3*pi/4),-R,R*cos(5*pi/4),0,R*cos(7*pi/4);...
         0,R*sin(pi/4),R,R*sin(3*pi/4),0,R*sin(5*pi/4),-R,R*sin(7*pi/4);...
         0,0,0,0,0,0,0,0];
r_UCA_halfwavelen = r_UCA/(lambda/2); % in half-wavelength

% Plot the UCA geometry
figure;
scatter3(r_UCA_halfwavelen(1,:),...
         r_UCA_halfwavelen(2,:),...
         r_UCA_halfwavelen(3,:),'MarkerFaceColor',[0,0.75,0.75]);
grid on;
xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
set(gca,'TickLabelInterpreter', 'latex');
xlabel('X-axis in half-wavelength','interpreter','latex');
ylabel('Y-axis in half-wavelength','interpreter','latex');
title('UCA Array Geometry','interpreter','latex');
set(gca,'FontSize',14);
view([0,90]);        
axis equal;

theta = 0:1:360; % azimuth
phi = -90:1:90;  % elevation

% desired beam direction in degrees
theta_steer = 90; % azimuth
phi_steer = 60;   % elevation
u_beam = [cos(deg2rad(theta_steer))*cos(deg2rad(phi_steer));...
          sin(deg2rad(theta_steer))*cos(deg2rad(phi_steer));...
          sin(deg2rad(phi_steer))];
w = exp(-j*2*pi/lambda*r_UCA'*u_beam);
phase_shifters = angle(conj(w));
phase_shifters = rad2deg(phase_shifters);
phase_shifters_deg = mod(phase_shifters,360);

% w_verify = exp(-1j*deg2rad(phase_shifters_deg));
% w = ones(N,1);

array_response = zeros(length(phi),length(theta));
X = zeros(1,length(theta)*length(phi));
Y = zeros(1,length(theta)*length(phi));
Z = zeros(1,length(theta)*length(phi));
Color = zeros(1,length(theta)*length(phi));
k = 0;
for i = 1:length(theta)
    for j = 1:length(phi)
        u = [cos(deg2rad(theta(i)))*cos(deg2rad(phi(j)));...
             sin(deg2rad(theta(i)))*cos(deg2rad(phi(j)));...
             sin(deg2rad(phi(j)))];  % directional vector
        % Calculate the gain of the array
        array_response(j,i) = abs(w'*exp(-1j*2*pi/lambda*r_UCA'*u));
%         array_response(j,i) = abs(sum(exp(1j*deg0rad(phase_shifters_deg)).*exp(-1j*2*pi/lambda*r_UCA'*u)));
        % Create 3D polar plot in x-y-z
        k = k + 1;
        X(1,k) = array_response(j,i)*cos(deg2rad(theta(i)))*cos(deg2rad(phi(j)));
        Y(1,k) = array_response(j,i)*sin(deg2rad(theta(i)))*cos(deg2rad(phi(j)));
        Z(1,k) = array_response(j,i)*sin(deg2rad(phi(j)));
        Color(1,k) = array_response(j,i);
        % Create 3D polar plot in azi-elevation
        azi(1,k) = theta(i);
        ele(1,k) = phi(j);
        Z_azi_ele(1,k) = array_response(j,i);
    end
end

% Linear
figure;
plot(theta,array_response((length(phi)+1)/2,:),'r','LineWidth',2);
grid on;
xlim([min(theta), max(theta)]);
% ylim([0,5]);
set(gca,'TickLabelInterpreter', 'latex');
xlabel('Azimuth Angle - degrees','interpreter','latex');
ylabel('Gain','interpreter','latex');
title('UCA Array Pattern','interpreter','latex');
set(gca,'FontSize',14);

% Polar on the Azimuth Plane
figure;
polarplot(deg2rad(theta),array_response((length(phi)+1)/2,:),'r','LineWidth',2);
set(gca,'TickLabelInterpreter', 'latex');
title('UCA Array Pattern (Polar)','interpreter','latex');
set(gca,'FontSize',14);

% 3D Cartesian
figure
scatter3(X,Y,Z,1,Color);
colormap(jet);
set(gca,'TickLabelInterpreter', 'latex');
xlabel('X','interpreter','latex');
ylabel('Y','interpreter','latex');
zlabel('Z','interpreter','latex');
title('UCA Array Pattern (3D-Cartesian)','interpreter','latex');
shading interp;
set(gca,'FontSize',14);
xlim([-8,8]);
ylim([-8 8]);
zlim([-8 8]);

xticks([-8:4:8]);
yticks([-8:4:8]);
zticks([-8:4:8]);
view(120,20);


% 3D spherical
figure;
scatter3(X,Y,Z,1,Color);
colormap(jet);
hold on;
%suppressing original x, y, axis
set(gca, 'XTick',[], 'ZTick',[], 'YTick',[], 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
xlim([-N-0.5,N+0.5]);ylim([-N-0.5,N+0.5]);zlim([-N-0.5,N+0.5]);


Azirange = (0:0.01:2*pi);   % azimuth angular range
angl = Azirange;
% Create circles for zero elevation
xa = N * cos(angl);
ya = N * sin(angl);
za = zeros(size(xa));
% Plot the polar axis
line(xa,ya,za,'Color','k','LineWidth',0.5);
for n = 1:1:N
    xa = (N-n) * cos(angl);
    ya = (N-n) * sin(angl);
    za = zeros(size(xa));
    % Plot the polar axis
    line(xa,ya,za,'Color','k','LineWidth',0.5); 
end
% Create tic marks at zero elevation
ts = 180/30;
ta = pi/ts * (round(0*ts/pi):1:round(2*pi*ts/pi));
tr = N  * [1; 1.03; 1.1];
xt = tr * cos(ta);
yt = tr * sin(ta);
zt = zeros(1,length(ta));
% Draw the tick marks
line(xt(1:2,:),yt(1:2,:),[zt; zt],'Color','k');
% Add tick labels
for k = 1:length(ta)-1
    text(xt(3,k),yt(3,k),zt(k),num2str(ta(k)*180/pi),'FontSize',12,'interpreter','latex');
end

Elerange = (0:0.01:2*pi);   % elevation angular range
angl = Elerange;
% Create circles for zero azimuth
xa = N * cos(angl);
ya = zeros(size(xa));
za = N * sin(angl);
line(xa,ya,za,'Color','k','LineWidth',0.5); 
% Create tic marks at zero azimuth
ts = 180/30;
ta = pi/ts * (round(0*ts/pi):1:round(pi*ts/pi)) - pi/2;
tr = N  * [1; 1.03; 1.1];
xt = tr * cos(ta);
yt = zeros(1,length(ta));
zt = tr * sin(ta);
% Draw the tick marks
line(xt(1:2,:),[yt;yt],zt(1:2,:),'Color','k');
% Add tick labels
for k = 1:length(ta)
    text(xt(3,k),yt(k),zt(3,k),num2str(ta(k)*180/pi),'FontSize',12,'interpreter','latex');
end

% Find x,y locations of radial labels
ta = 45*pi/180;
tr = (0:1:N);
xt = tr * cos(ta);
yt = tr * sin(ta);
zt = zeros(1,length(tr));
% Add radial labels to the plot
for k = 1:length(tr)
    text(xt(k),yt(k),zt(k),num2str(tr(k)),'FontSize',12,'interpreter','latex');
end
line(xt,yt,zt,'Color','k','LineWidth',1);
view(120,20);
title('UCA Array Pattern (3D-Spherical)','interpreter','latex');
set(gca,'FontSize',14);
axis square