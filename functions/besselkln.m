function l = besselkln(nu,z)
% Log modified Bessel of second kind
l = log(besselk(nu,z));
if any(l==-Inf)
    if abs(nu)>50
        z0 = min(z(l==-Inf));
        z0 = find(log(besselk(nu,1:z0))>-Inf,1,'last');
        l(l==-Inf) = z0-z(l==-Inf)-0.5*log(z(l==-Inf)/z0)+log(besselk(nu,z0));
    else
        l(l==-Inf) = -z(l==-Inf)-0.5*log(z(l==-Inf))+0.5*log(pi/2);
    end
end
    
if any(l==Inf)
    if nu~=0
        l(l==Inf) = gammaln(abs(nu))+abs(nu)*(log(2)-log(z(l==Inf)))-log(2);
    else
        l(l==Inf) = log(-log(z(l==Inf)/2)-eulergamma);
    end
end
end

% 
% function k = approxlb(nu,z)
% k = -z-0.5*log(z)+0.5*log(pi/2)+log(1+(4*nu^2-1)*z.^(-1)/8+...
%     (4*nu^2-1)*(4*nu^2-9)*z.^(-2)/128+...
%     (4*nu^2-1)*(4*nu^2-9)*(4*nu^2-25)*z.^(-3)/3072);