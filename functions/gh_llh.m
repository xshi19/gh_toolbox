function [llh,lambda_,eta_,delta_,chi_,psi_] = ...
    gh_llh(sample,mu,gamma,Sigma,lambda,chi,psi)
% Loglikelihood of GH
[mu,gamma,Sigma,lambda,chi,psi] = ...
    gh_reg(mu,gamma,Sigma,lambda,chi,psi,4);
[n,d] = size(sample);
R = chol(Sigma);
sample_ = (sample-repmat(mu',n,1))/R;
gamma_ = R'\gamma;
lambda_ = lambda-d/2; %scalar
psi_ = psi+gamma_'*gamma_;
chi_ = chi+sum(sample_.^2,2);
delta_ = sqrt(chi_/psi_);
eta_ = sqrt(chi_*psi_);

if n>100
    lbes = besselkln2(lambda_,eta_);
else
    lbes = besselkln(lambda_,eta_);
end
    
llh = lambda/2*(log(psi/chi))...
    - log(besselk(lambda,sqrt(chi*psi)))...
    + lbes + sample_*gamma_...
    + lambda_/2*(log(chi_/psi_));
end