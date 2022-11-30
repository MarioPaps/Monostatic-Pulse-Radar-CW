function pdf = fSwerling34pdf(sigma,x,type)
%FSWERLING34PDF 
%   Inputs:   sigma: mean value of RCS in m^2
%             x: RCS or Amplitude
%             type: the RCS or Amplitude model
%   Outputs:  pdf: is the pdf at values of x
switch type
    case 'RCS'
%         pdf = 4*x/sigma^2.*exp(-2*x/sigma);           
        pdf = 4/sigma*chi2pdf((4/sigma)*x,4);           % chi-sq DOF=4
    case 'Amplitude'
%         pdf = 8*x.^3/sigma^2.*exp(-2*x.^2/sigma);      
        pdf = 4/sigma^4*x.^2.*raylpdf(2/sqrt(sigma)*x,1); 
end
end

