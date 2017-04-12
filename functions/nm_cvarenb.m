function [N,p] = nm_cvarenb(w,alpha,mu,gamma,Sigma,y,z)

[~,d,H] = nm_portcvar(w,alpha,gamma,Sigma,y,z);
d = d-mu;
H2 = 2*(d*d') + 2*(w'*d)*H;
[N, p] = enbCov(w, H2);

end