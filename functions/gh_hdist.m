function d = gh_hdist(mu1,gamma1,Sigma1,lambda1,chi1,psi1,...
    mu2,gamma2,Sigma2,lambda2,chi2,psi2)
% Hellinger distance between two GH (joint) distributions

dm = mu1-mu2;
dg = gamma1-gamma2;
Sigma_ = (Sigma1+Sigma2)/2;

lambda_ = (lambda1+lambda2)/2;
chi_ = (chi1+chi2)/2+dm'*(Sigma_\dm);
psi_ = (psi1+psi2)/2+dg'*(Sigma_\dg);

c1 = (psi1/chi1)^(lambda1/2)/besselk(lambda1,sqrt(chi1*psi1));
c2 = (psi2/chi2)^(lambda2/2)/besselk(lambda2,sqrt(chi2*psi2));
c_ = (psi_/chi_)^(lambda_/2)/besselk(lambda_,sqrt(chi_*psi_));

d = 1-det(Sigma1*Sigma2)^(1/4)/det(Sigma_)^(1/2)*...
    sqrt(c1*c2)/c_*exp(-dm'*(Sigma_\dg)/4);

end
    