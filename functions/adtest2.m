function [H,adstat,critvalue] = adtest2 (sample1, sample2, alpha)
% Sample1 and sample2 are the samples to be compared. They must be vectors of a size greater than 8. Alpha specifies the
% allowed error. If alpha is not specified, a default value of
% 0.05 for alpha is used. Alpha must be either 0.01, 0.05 or 0.1.
% RETURN: H gives the statistical decision. H = 0: samples are not significantly different. H = 1: sample1 and sample2 are significantly different (i.e. do not arise from the same underlying distribution).
% adstat returns the ADstatistic of the comparison of the two samples. If adstat is greater than the critical value,
% the two samples are significantly different.
% critvalue returns the critical value for the alpha used (c) Sonja Engmann 2007

if nargin < 2
    error('Call adtest2 with at least two input arguments'); 
end
if nargin < 3
    alpha = 0.05;
end
% Assignment of critical value depending on alpha 
if alpha == 0.01
    critvalue = 3.857;
elseif alpha == 0.05
    critvalue = 2.492;
elseif alpha == 0.1
    critvalue = 1.933;
else
    error('Alpha must be either 0.01, 0.05 or 0.1.');
end

samplecomb = sort([sample1(:)' sample2(:)']); 


m = interp1(sort(sample1(:)),(1:length(sample1))',samplecomb(1:end-1),'previous','extrap');
m(isnan(m)) = 0;
n = 1:length(samplecomb)-1;
ad = sum((m*length(samplecomb)-n*length(sample1)).^2./(n.*(length(samplecomb)-n)));


adstat = ad/(length(sample1)*length(sample2)); 
if adstat > critvalue
    H = 1; 
else
    H = 0;
end