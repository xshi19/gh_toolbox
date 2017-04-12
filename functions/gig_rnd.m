function z = gig_rnd(n,lambda,chi,psi,lb,ub,dx,u,method)
% generate GIG random numbers    
if nargin<=4
    m = gig_moment(1,lambda,chi,psi);
    v = gig_moment(2,lambda,chi,psi)-m^2;
    lb = 1e-10;
    ub = m+20*sqrt(v);
    dx = (ub-lb)/(2*n+1);
    u = unifrnd(0,1,n,1);
    method = 'pchip';
end
x = (lb:dx:ub)';
p = gig_pdf(x,lambda,chi,psi);
p(p<=0) = 0;
y = cumsum(p)/sum(p);
[~,ind] = unique(y);
z = interp1(y(ind),x(ind),u,method);
end