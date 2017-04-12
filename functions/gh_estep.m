function [a,b,c] = gh_estep(lambda_,eta_,delta_)
% E-step of the EM algorithm
% computes conditional expectations of Y^(-1), Y and log(Y).
a = besselkratio2(lambda_-1,lambda_,eta_);
b = besselkratio2(lambda_+1,lambda_,eta_);
c = (besselkratio2(lambda_+1e-10,lambda_,eta_)-...
    besselkratio2(lambda_-1e-10,lambda_,eta_))/2e-10;
% aux0 = besselkln(lambda_,eta_);
% aux1 = besselkln(lambda_-1,eta_);
% aux2 = besselkln(lambda_+1,eta_);
% aux3 = besselkln(lambda_+1e-10,eta_);
% aux4 = besselkln(lambda_-1e-10,eta_);
% 
% a = exp(aux1-aux0);
% b = exp(aux2-aux0);
% c = (exp(aux3-aux0)-exp(aux4-aux0))/2e-10;

% c = (aux3-aux4)./aux0/2e-10;
% 
% la = logical((a>0).*(a<Inf));
% lb = logical((b>0).*(b<Inf));
% lc = logical((c>-Inf).*(c<Inf));
% 
% a(~la) = 0.5*(sqrt(((2*lambda_-1)./eta_(~la)).^2+4)-...
%     (2*lambda_-1)./eta_(~la));
% b(~lb) = 2./(sqrt(((2*lambda_+1)./eta_(~lb)).^2+4)-...
%     (2*lambda_+1)./eta_(~lb));
% 
% if sum(lc)~=0
%     c(~lc) = interp1(eta_(lc),c(lc),eta_(~lc),'previous','extrap');
%     lc = logical((c>-Inf).*(c<Inf));
%     c(~lc) = interp1(eta_(lc),c(lc),eta_(~lc),'pchip','extrap');
% end

a = a./delta_;
b = b.*delta_;
c = c + log(delta_);

end