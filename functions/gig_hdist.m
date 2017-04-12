function d = gig_hdist(lambda1,chi1,psi1,lambda2,chi2,psi2)
% Hellinger distance of GIG
lambda_ = (lambda1+lambda2)/2;
chi_ = (chi1+chi2)/2;
psi_ = (psi1+psi2)/2;

c1 = (psi1/chi1)^(lambda1/2)/besselk(lambda1,sqrt(chi1*psi1));
c2 = (psi2/chi2)^(lambda2/2)/besselk(lambda2,sqrt(chi2*psi2));
c_ = (psi_/chi_)^(lambda_/2)/besselk(lambda_,sqrt(chi_*psi_));
d = 1 - sqrt(c1*c2)/c_;

end
   