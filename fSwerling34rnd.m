function rnd = fSwerling34rnd(sigma,n,type)
%FSWERLING34RND
%   Inputs:   sigma: mean value of RCS in m^2
%             n: number of samples
%             type: the RCS or Amplitude model
%   Outputs:  rnd: the random number following specific model (pdf)
switch type
    case 'RCS'
        rnd = sigma/4*chi2rnd(4,[1,n]);
    case 'Amplitude'
        rnd = sqrt(sigma/4*chi2rnd(4,[1,n]));
end
end

