function [mu,gamma,Sigma,lambda,chi,psi] = ...
    gh_reg(mu,gamma,Sigma,lambda,chi,psi,reg)
% Regulation of GH parameters

if reg == 1
    % chi = 1
    gamma = gamma*chi;
    Sigma = Sigma*chi;
    psi = psi*chi;
    chi = 1;
elseif reg == 2
    % psi = 1
    gamma = gamma/psi;
    Sigma = Sigma/psi;
    chi = chi*psi;
    psi = 1;
elseif reg == 3
    % chi = psi
    gamma = gamma*sqrt(chi/psi);
    Sigma = Sigma*sqrt(chi/psi);
    chi = sqrt(chi*psi);
    psi = chi;
elseif reg == 4
    % det(Sigma) = 1
    if det(Sigma)==0
        c = eps^(25/size(Sigma,1));
        d = det(Sigma/c)^(1/size(Sigma,1));
        d = d*c;
    elseif det(Sigma)==Inf
        c = eps^(25/size(Sigma,1));
        d = det(Sigma*c)^(1/size(Sigma,1));
        d = d/c;
    else
        d = det(Sigma)^(1/size(Sigma,1));
    end
    gamma = gamma/d;
    Sigma = Sigma/d;
    chi = chi*d;
    psi = psi/d;
end

end