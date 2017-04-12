function l = besselkln2(nu,z)
% Log modified Bessel of second kind for a sequence z
l = log(besselk(nu,z));
if any(l==-Inf)
    z0 = max(z(l>-Inf));
    l(l==-Inf) = z0-z(l==-Inf)-0.5*log(z(l==-Inf)/z0)+log(besselk(nu,z0));
end
    
if any(l==Inf)
    z0 = min(z(l<Inf));
    if isempty(z0)
        l(l==Inf) = -abs(nu)*log(z(l==Inf))+(nu-1)*log(2)+gammaln(abs(nu));
    else
        l(l==Inf) = -abs(nu)*log(z(l==Inf)/z0)+log(besselk(nu,z0));
    end
end
end