function [mu,gamma,Sigma,lambda,chi,psi,llh,F,D] = ...
    gh_fac(sample,mu0,gamma0,F0,D0,lambda0,chi0,psi0,maxiter,distr,fixt,thres)
% Factor analysis for Generalized Hyperbolic distributions
%
% distr: distribution type. Avalible distributions are
%   'GH': Generalized Hyperbolic distribution
%   'GHl': Generalized Hyperbolic distribution with fixed lambda
%   'NIG': Normal inverse Gaussian distribution
%   'VG': Variance gamma distribtuion
%   'T': Skewed T distribution
% fixt: 
%   0: Do not fix tail parameters (lambda,chi,psi)
%   1: Fix tail parameters (lambda,chi,psi)

[n,dim] = size(sample);
S4 = mean(sample)';

if nargin<12
    thres = 1e-5;
end
if nargin<11
    fixt = 0;
end
if nargin<10
    distr = 'GH';
end
if nargin<9
    maxiter = 100;
end
if nargin<2
    r = floor(dim/2);
    mu0 = mean(sample)';
    gamma0 = zeros(dim,1);
    lambda0 = -0.5;
    chi0 = 1;
    psi0 = 1;
    Sigma0 = cov(sample)/gig_moment(1,lambda0,chi0,psi0);
    [V,L] = eigs(Sigma0,r);
    F0 = V*sqrt(L);
    D0 = diag(Sigma0-F0*F0');
end

r = size(F0,2);
mu = mu0;
gamma = gamma0;
F = F0;
D = D0;
Sigma = F0*F0'+diag(D0);
lambda = lambda0;
chi = chi0;
psi = psi0;

options = optimoptions('fmincon','GradObj','on',...
    'Display','off','Algorithm','interior-point');

switch distr
    case 'GH'
        lb = [-20;0;0];
        ub = [20;Inf;Inf];
        Aeq = [];
        beq = [];
    case 'GHl'
        lb = [-20;0;0];
        ub = [20;Inf;Inf];
        Aeq = [1,0,0];
        beq = lambda0;
    case 'NIG'
        lb = [-20;0;0];
        ub = [20;Inf;Inf];
        lambda = -0.5;
        Aeq = [1,0,0];
        beq = -0.5;
    case 'VG'
        chi = 0;
        lambda = max(1e-5,lambda);
        lb = [0;0;0];
        ub = [20;Inf;Inf];
        Aeq = [0,1,0];
        beq = 0;
end

llh = zeros(maxiter,1);

% set det(Sigma)=1
[~,d] = setdet1(Sigma);
gamma = gamma/d;
Sigma = Sigma/d;
chi = chi*d;
psi = psi/d;
F = F/sqrt(d);
D = D/d;

for iter=1:maxiter 
    % iter
    % Main iteration
    [llh_,lambda_,eta_,delta_] = ...
        gh_llh(sample,mu,gamma,Sigma,lambda,chi,psi);
    llh(iter) = mean(llh_);
    
    if iter>1
        if abs(llh(iter)-llh(iter-1))<thres
            break;
        end
    end
    
    % a, b, c: n times 1
    [a,b,c] = gh_estep(lambda_,eta_,delta_);
    
    % sufficient statistics
    beta = Sigma\F;
    S1 = mean(a);
    S2 = mean(b);
    S3 = mean(c);
    S5 = sample'*a/n;
    S6 = repmat(sqrt(a),1,dim).*sample;
    S6 = S6'*S6/n;
    S7 = (S6-S5*mu'-S4*gamma')*beta;
    S8 = beta'*(S5-mu*S1-gamma);
    S9 = beta'*(S4-mu-gamma*S2);
    S10 = eye(r)-beta'*F+beta'*(S6-S5*mu'-mu*S5'+mu*mu'*S1-...
        (S4-mu)*gamma'-gamma*(S4-mu)'+gamma*gamma'*S2)*beta;
    
    % M step
    u = S10\S8;
    v = S10\S9;
    T1 = S8'*u-S1;
    T2 = S9'*u-1;
    T3 = S9'*v-S2;
    T4 = S7*u-S5;
    T5 = S7*v-S4;
    
    mu = (T2*T5-T3*T4)/(T2^2-T1*T3);
    gamma = (T2*T4-T1*T5)/(T2^2-T1*T3);
    F = (S7-mu*S8'-gamma*S9')/S10;
    D = diag(S1*(mu*mu')+S2*(gamma*gamma')-S4*gamma'-gamma*S4'-S5*mu'-mu*S5'+...
        S6-S7*F'-F*S7'+F*S8*mu'+mu*(F*S8)'+F*S9*gamma'+gamma*(F*S9)'+...
        F*S10*F'+mu*gamma'+gamma*mu');
    if fixt==0
        try
            [lambda,chi,psi] = ...
                gig_mle(S1,S2,S3,[lambda;chi;psi],Aeq,beq,lb,ub,options);
        catch
            disp('Error in optimization.')
        end
    end
    
    Sigma = F*F'+diag(D);
    [Sigma,d] = setdet1(Sigma);
    gamma = gamma/d;
    chi = chi*d;
    psi = psi/d;
    F = F/sqrt(d);
    D = D/d;
end

llh(iter+1:end) = llh(iter);

end

