function p = gig_pdf(x,lambda,chi,psi)
% pdf of Generalized Inverse Gaussian (GIG)
c = (psi/chi)^(lambda/2)/2/besselk(lambda,sqrt(chi*psi));
if c==Inf
    c = (psi/2)^lambda/gamma(lambda);
end

p = c*x.^(lambda-1).*exp(-0.5*(chi*x.^(-1)+psi*x)).*(x>0);
end