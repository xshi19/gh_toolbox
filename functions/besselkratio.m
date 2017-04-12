function r = besselkratio(nu1,nu2,z)
b1 = besselk(nu1,z);
b2 = besselk(nu2,z);
r = b1./b2;

loc = ~logical((b1<Inf).*(b2<Inf));
if any(loc)
    r(loc) = exp(gammaln(abs(nu1))-gammaln(abs(nu2)))*(z(loc)/2).^(abs(nu2)-abs(nu1));
end

if any(b2==0)
    z0 = min(z(b2==0));
    z0 = find(besselk(nu1,1:z0)>0,1,'last');
    r(b2==0) = besselk(nu1,z0)./besselk(nu2,z0);
end
