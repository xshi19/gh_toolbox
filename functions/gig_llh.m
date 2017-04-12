function [l,g] = gig_llh(lambda,chi,psi,S1,S2,S3)
% -log likelihood of GIG
eta = sqrt(chi.*psi);
delta = sqrt(chi./psi);
l = 0.5*(chi*S1+psi*S2)-(lambda-1)*S3+...
    lambda*log(delta)+log(besselk(lambda,eta));

if nargout > 1
%     aux0 = besselkln(lambda,eta);
%     aux1 = besselkln(lambda-1,eta);
%     aux2 = besselkln(lambda+1,eta);
%     aux3 = besselkln(lambda+1e-10,eta);
%     aux4 = besselkln(lambda-1e-10,eta);
%     g = [-S3+log(delta)+(exp(aux3-aux0)-exp(aux4-aux0))/2e-10;...
%         0.5*(S1-exp(aux1-aux0)./delta);
%         0.5*(S2-exp(aux2-aux0).*delta)];
    g = [-S3+log(delta)+(besselkratio(lambda+1e-10,lambda,eta)-...
        besselkratio(lambda-1e-10,lambda,eta))/2e-10;...
        0.5*(S1-besselkratio(lambda-1,lambda,eta)./delta);
        0.5*(S2-besselkratio(lambda+1,lambda,eta).*delta)];
end
end