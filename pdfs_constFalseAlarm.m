clc,clear,close all;
%% Written by Yunhao Liu
%% Last modified: 02/03/2021
%% Constant False Alarm
%% Find the value of threshold with a pre-defined probability of false alarm
Vol_max = 40;
Vol_step = 0.01;
Voltage = (0:Vol_step:Vol_max);  % voltage range
datapoints = length(Voltage);

ProbFA = 0.9; % Pre-defined Probability of False Alarm

% Note: quadrature noise representation + Envelope detector Rx
Noise_var = 5;  % noise power
Noise_mean = 0; % mean
Noise = raylpdf(Voltage,sqrt(Noise_var/2));


% Note: clutter -- Weibull distribution
% Weibull shape parameter (Clutter_shape = 1: exponential; = 2: Rayleigh)
Clutter_shape = 1.2; 
Clutter_var = 1; % clutter power
Clutter_scale = sqrt(Clutter_var/(gamma(1+2/Clutter_shape)-(gamma(1+1/Clutter_shape))^2)); % clutter Weibull scale parameter
Clutter = wblpdf(Voltage,Clutter_scale,Clutter_shape);  % clutter pdf
Clutter_mean = Clutter_scale*gamma(1+1/Clutter_shape);


% Control of SNR with fixed noise and clutter power
SNR = 10; % Signal-to-Noise Ratio in dB

Noise_var_dB = 10*log10(Noise_var);
Signal_var_dB = SNR + Noise_var_dB;
Signal_var = 10^(Signal_var_dB/10);

% Signal-to-Noise-plus-Clutter Ratio in dB
SNCR = 10*log10(Signal_var/(Noise_var+Clutter_var));

% Signal/Target -- Rayleigh
Signal_mean = 0;
Signal = raylpdf(Voltage,sqrt(Signal_var/2));     % signal pdf


% Plot of three pdfs
figure;
plot(Voltage,Noise,'Color','k','LineWidth',2);
hold on;
plot(Voltage,Clutter,'Color','m','LineWidth',2);
plot(Voltage,Signal,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
set(gca,'TickLabelInterpreter', 'latex');
legend('Noise pdf (Rayleigh)',...
       'Clutter pdf (Weibull)',...
       'Signal pdf (Rayleigh)',...
       'interpreter','latex');
ylim([0 0.6]);
xlim([0 18]);
xlabel('Amplitude','interpreter','latex');
ylabel('pdf','interpreter','latex');
set(gca,'FontSize',14);


% Noise plus Clutter --> convolution of two pdfs
Noise_Clutter = conv(Noise,Clutter)*Vol_step;
% % convolution using fft and ifft
% Noise_ext = [Noise zeros(1,datapoints-1)]; % zero padding
% Clutter_ext = [Clutter zeros(1,datapoints-1)]; % zero padding
% Noise_Clutter = ifft(fft(Noise_ext).*fft(Clutter_ext))*Vol_step;


% Noise plus Clutter plus Signal
Noise_Clutter_Signal = conv(Noise_Clutter,Signal)*Vol_step;
% % convolution using fft and ifft
% NC_ext = [Noise_Clutter, zeros(1,datapoints-1)];
% Signal_ext = [Signal, zeros(1,length(Noise_Clutter)-1)];
% Noise_Clutter_Signal = ifft(fft(NC_ext).*fft(Signal_ext))*Vol_step

%cumtrapz(noise_clutter) computes the cdf of the noise+clutter pdf
Vol_idx = find(cumtrapz(Noise_Clutter)*Vol_step<(1-ProbFA));
thres_FA = Voltage(max(Vol_idx)); % Find the corresponding threshold
% Verification
ProbFA_verification = sum(Noise_Clutter(1,length(0:Vol_step:thres_FA):end))*Vol_step;


% Plot of two likelihood functions, i.e. with and without target present
figure;
vol_temp0 = (thres_FA:Vol_step:2*Vol_max);
curve1 = Noise_Clutter(1,length(0:Vol_step:thres_FA):end);
curve0 = zeros(1,length(vol_temp0));
x_temp_plot = [vol_temp0,fliplr(vol_temp0)];
fill(x_temp_plot,[curve1,fliplr(curve0)],'y','facealpha',1);
hold on;grid on;
vol_temp0 = (0:Vol_step:thres_FA);
curve1 = Noise_Clutter_Signal(1,1:length(vol_temp0));
curve0 = zeros(1,length(vol_temp0));
x_temp_plot = [vol_temp0,fliplr(vol_temp0)];
fill(x_temp_plot,[curve1,fliplr(curve0)],'g','facealpha',1);
plot((0:Vol_step:2*Vol_max),Noise_Clutter,'LineWidth',2);
hold on;
plot((0:Vol_step:3*Vol_max),Noise_Clutter_Signal,'r','LineWidth',2);
xline(thres_FA,'-',{'Threshold'},'LineWidth',2);% line for detection threshold
ylim([0 0.3]);
xlim([0 25]);
set(gca,'TickLabelInterpreter', 'latex');
xlabel('Amplitude','interpreter','latex');
ylabel('pdf','interpreter','latex');
legend(['Probability of False Alarm',' ', num2str( ProbFA)],...
       'Probability of a Miss',...
       'Noise pdf * Clutter pdf',...
       'Noise pdf * Clutter pdf * Signal pdf',...
       'interpreter','latex');
set(gca,'FontSize',14);