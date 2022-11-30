

function [threshold] = specifyThreshold(noise,tasknum)

    ProbFA=0.001;
    noise_real= abs(noise);
    Pn_est= (1/length(noise))* (noise)* ctranspose(noise);
    vol_max= max(noise_real);
    vol_min= min(noise_real);
    [N,edges]= histcounts(noise_real,'Normalization','pdf');
    vol_step= edges(2)-edges(1);
    if(tasknum~=7)
        voltage_range= vol_min:vol_step:vol_max;
    else
        voltage_range= linspace(vol_min,vol_max,length(noise_real));
    end

%     numbins= fix((vol_max*2-vol_min)/2e-2);
    noise_pdf= raylpdf(noise_real,sqrt(Pn_est/2));
%     vol_idx = find(cumtrapz(noise_pdf)*vol_step<(1-ProbFA));
%     threshold= voltage_range(max(vol_idx));
    vol_idx= find(cumtrapz(noise_pdf)*vol_step<(1-ProbFA),1,'last');
    threshold= voltage_range(vol_idx);


%     ProbFA_verification = sum(noise_pdf(1,length(0:vol_step:threshold):end))*vol_step



%     %compute the cdf of the noise distribution
%    
%     threshold=abs( raylinv(1-P_fa,sqrt(Pn_est/2)) );
% 
% %     cdfnoise=raylcdf(abs(noisysig),sqrt(Pn_est/2));
% %     [cdfnoise_ord,idx]= sort(cdfnoise,'ascend');
% %     value= cdfnoise_ord(length(cdfnoise_ord)-13);
% %     pos= find(cdfnoise_ord==value);
% %     original_pos= idx(pos);
% %     threshold= noisysig(original_pos);
% 
%     
% %     threshold=0.1;
% %     pdf_n= h.Values;
%     ds= abs(edges(2)-edges(1));
%     v=0:ds: max(noisysig);
%     cdf=0;
%     for ind=1:length(pdf_n)
%         cdf=cdf+ ds* pdf_n(ind);
%           if(cdf< (1-P_fa))
%                 threshold= 0.1;
%           end
%           endc


end