function rnd = fSwerling12rnd(sigma,n,type)
%FSWERLING12RND
%   Inputs:   sigma: mean value of RCS in m^2
%             n: number of samples
%             type: the RCS or Amplitude model
%   Outputs:  rnd: the random number following specific model (pdf)
switch type
    case 'RCS'
        rnd = exprnd(sigma,[1,n]);
    case 'Amplitude'
%         rnd = sqrt(exprnd(sigma,[1,n]));
        rnd = raylrnd(sqrt(sigma/2),[1,n]);
end
end

