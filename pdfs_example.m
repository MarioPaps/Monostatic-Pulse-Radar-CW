clc,clear,close all;

L = 10000;
Pn = 0.8;
sigma = sqrt(Pn/2);

% Generate random numbers
option = 'a';
switch option
    case 'a'
        nreal = randn(1,L);
        nimag = randn(1,L);
        nsamples = sqrt(Pn/2)*(nreal+1j*nimag);
    case 'b'
        nreal = normrnd(0,sigma,1,L);
        nimag = normrnd(0,sigma,1,L);
        nsamples = nreal+1j*nimag;
end

% Gaussian
% estimate pdf
nreal = real(nsamples);
NBins = fix((max(nreal)*2-min(nreal))/2e-2);
[Nreal_pdf_hat,edges] = histcounts(nreal,NBins,'Normalization','pdf');
edges = mean([edges(1:end-1);edges(2:end)]);
% theoretical pdf
Volt = (-3:0.01:3);
Nreal_pdf = normpdf(Volt,0,sigma);

figure;
plot(edges,Nreal_pdf_hat,'LineWidth',1);
hold on; grid on;
plot(Volt,Nreal_pdf,'Color','k','LineWidth',2);
title('Gaussian distribution');
xlabel('Amplitude');
ylabel('pdf');
set(gca,'FontSize',12);


% Rayleigh
% estimate pdf
Nmag = abs(nsamples);
NBins = fix((max(Nmag)*2-min(Nmag))/2e-2);
[N_pdf_hat,edges] = histcounts(Nmag,NBins,'Normalization','pdf');
edges = mean([edges(1:end-1);edges(2:end)]);
% theoretical pdf
Mag = (0:0.01:3);
N_pdf = raylpdf(Mag,sqrt(Pn/2));

figure;
plot(edges,N_pdf_hat,'LineWidth',1);
hold on; grid on;
plot(Mag,N_pdf,'Color','k','LineWidth',2);
title('Rayleigh distribution');
xlabel('Magnitude');
ylabel('pdf');
set(gca,'FontSize',12);