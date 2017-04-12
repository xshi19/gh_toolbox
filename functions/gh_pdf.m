function p = gh_pdf(x,mu,gamma,Sigma,lambda,chi,psi)
% pdf of Generalized Hyperbolic (GH) distribution
[n,d] = size(x);

% [mu,gamma,Sigma,lambda,chi,psi] = gh_reg(mu,gamma,Sigma,lambda,chi,psi,4);
R = chol(Sigma);
x_ = (x-repmat(mu',n,1))/R;
gamma_ = R'\gamma;
lambda_ = lambda-d/2;
chi_ = chi+sum(x_.^2,2);
psi_ = psi+gamma_'*gamma_;

c = (psi/chi)^(lambda/2)/(2*pi)^(d/2)/...
    sqrt(det(Sigma))/besselk(lambda,sqrt(chi*psi));
p = besselk(lambda_,sqrt(chi_*psi_)).*...
    exp(x_*gamma_).*sqrt(psi_./chi_).^(d/2-lambda);
p = c.*p;

end