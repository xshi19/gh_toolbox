function [S1,S2,S3,S4,S5,S6] = param2sufstat(mu,gamma,Sigma,lambda,chi,psi)

S1 = gig_moment(-1,lambda,chi,psi);
S2 = gig_moment(1,lambda,chi,psi);
S3 = (gig_moment(1e-10,lambda,chi,psi) - ...
    gig_moment(-1e-10,lambda,chi,psi)) / 2e-10;
S4 = mu+gamma*S2;
S5 = mu*S1+gamma;
S6 = Sigma+S1*(mu*mu')+S2*(gamma*gamma')+mu*gamma'+gamma*mu';