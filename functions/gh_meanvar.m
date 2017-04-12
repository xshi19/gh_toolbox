function [m,C] = gh_meanvar(mu,gamma,Sigma,lambda,chi,psi)
a = gig_moment(1,lambda,chi,psi);
b = gig_moment(2,lambda,chi,psi)-a^2;
m = mu+a*gamma;
C = a*Sigma+b*(gamma*gamma');