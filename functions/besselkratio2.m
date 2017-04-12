function r = besselkratio2(nu1,nu2,z)
% Ratio of two besselk with different nu for a sequence z
inv = 0;
if abs(nu1)<abs(nu2)
    nu_ = nu1;
    nu1 = nu2;
    nu2 = nu_;
    inv = 1;
end

b1 = besselk(nu1,z);
b2 = besselk(nu2,z);
r = b1./b2;

loc = logical((b1<Inf).*(b2<Inf));
if any(~loc)
    z0 = min(z(loc));
    if isempty(z0)
        r(~loc) = exp(gammaln(abs(nu1))-gammaln(abs(nu2)))*...
            (z(~loc)/2).^(abs(nu2)-abs(nu1));
    else
        r(~loc) = exp(gammaln(abs(nu1))-gammaln(abs(nu2)))*...
            ((z(~loc)/2).^(abs(nu2)-abs(nu1))-(z0/2).^(abs(nu2)-abs(nu1)))+...
            besselk(nu1,z0)/besselk(nu2,z0);
    end
end

if any(b2==0)
    z0 = max(z(b2>0));
    r(b2==0) = besselk(nu1,z0)./besselk(nu2,z0);
end

if inv
    r = r.^(-1);
end

