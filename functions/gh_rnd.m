function [x,y] = gh_rnd(n,mu,gamma,Sigma,lambda,chi,psi)
d = length(mu);
y = gig_rnd(n,lambda,chi,psi);
% y = randraw('gig',[lambda,chi,psi],n);
x = repmat(mu',n,1)+repmat(y,1,d).*repmat(gamma',n,1)+...
    repmat(sqrt(y),1,d).*mvnrnd(zeros(1,d),Sigma,n);
end

