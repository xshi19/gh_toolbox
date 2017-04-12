function [v,d,H] = nm_portcvar(w,alpha,gamma,Sigma,y,z)

sigma = sqrt(w'*Sigma*w);
[v,dg,ds,dgg,dgs,dss] = nm_cvar1(alpha,w'*gamma,sigma,y,z);

d = dg*gamma+ds*Sigma*w/sigma;
H = gamma*gamma'*dgg+(gamma*(w'*Sigma)+(Sigma*w)*gamma')*dgs/sigma...
    +(Sigma*w)*(Sigma*w)'*dss/sigma^2+(Sigma/sigma-...
    (Sigma*w)*(Sigma*w)'/sigma^3)*ds;

end
