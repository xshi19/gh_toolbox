function m = gig_moment(n,lambda,chi,psi)
% compute the n-th moment of GIG
% m = exp(log(chi/psi)*n/2).*...
%     besselk(lambda+n,sqrt(chi*psi))./besselk(lambda,sqrt(chi*psi));
delta = sqrt(chi/psi);
eta = sqrt(chi*psi);
m = delta.^n.*exp(besselkln(lambda+n,eta)-besselkln(lambda,eta));
end