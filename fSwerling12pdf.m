function pdf = fSwerling12pdf(sigma,x,type)
%FSWERLING12PDF 
%   Inputs:   sigma: mean value of RCS in m^2
%             x: RCS or Amplitude
%             type: the RCS or Amplitude model
%   Outputs:  pdf: is the pdf at values of x
switch type
    case 'RCS'
%         pdf = 1/sigma*exp(-x/sigma);                  
        pdf = exppdf(x,sigma);                        % exponential DOF=4
%         pdf = 2/sigma*chi2pdf((2/sigma)*x,2);         % chi-sq DOF=2        
    case 'Amplitude'
%         pdf = 2*x/sigma.*exp(-x.^2/sigma);            
        pdf = raylpdf(x,sqrt(sigma/2));               % Rayleigh
end
end

