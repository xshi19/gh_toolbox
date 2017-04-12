function [cv,dg,ds,dgg,dgs,dss] = nm_cvar1(alpha,gamma,sigma,y,z)
% CVaR of 1 dimensional normal mixture
n = floor(length(z)*alpha);
xs = sort(gamma*y+sigma*sqrt(y).*z);
cv = -mean(xs(1:n));
x_ = (xs(n)-gamma*y)./(sigma*sqrt(y));

dg = -mean(y.*normcdf(x_))/alpha;
ds = (cv-gamma*dg)/sigma;

dg_ = -mean(sqrt(y).*normpdf(x_))/mean(normpdf(x_)./sqrt(y)); % dVaR/dgamma

dgg = mean(sqrt(y).*normpdf(x_).*(dg_+y))/alpha/sigma;
dgs = -gamma*dgg/sigma;
dss = -gamma*dgs/sigma;

end