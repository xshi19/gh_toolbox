function [f,g] = nm_cvarobj(w,w0,alpha,m,l,c,gamma,Sigma,y,z)
[v,d] = nm_portcvar(w,alpha,gamma,Sigma,y,z);
f = -w'*m+l*v+c*norm(w-w0,1);
g = -m+l*d+c*sign(w-w0);
end