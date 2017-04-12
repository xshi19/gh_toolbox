function [mu,gamma,Sigma,lambda,chi,psi,llh] = ...
    gh_mcecm(sample,mu0,gamma0,Sigma0,lambda0,chi0,psi0,maxiter,distr,fixt,thres)
% MCECM algorithm for Generalized Hyperbolic distributions
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

if nargin<11
    thres = 1e-5;
end
if nargin<10
    fixt = 0;
end
if nargin<9
    distr = 'GH';
end
if nargin<8
    maxiter = 100;
end
if nargin<2
    mu0 = mean(sample)';
    gamma0 = zeros(dim,1);
    lambda0 = -0.5;
    chi0 = 1;
    psi0 = 1;
    Sigma0 = cov(sample)/gig_moment(1,lambda0,chi0,psi0);
end

mu = mu0;
gamma = gamma0;
Sigma = Sigma0;
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
    case 'T'
        lambda = -chi/2;
        lb = [-20;0;0];
        ub = [20;Inf;Inf];
        Aeq = [2,1,0];
        beq = 0;
end

llh = zeros(maxiter,1); % log likelihood
[mu,gamma,Sigma,lambda,chi,psi] = ...
    gh_reg(mu,gamma,Sigma,lambda,chi,psi,4);

for iter=1:maxiter 
    % Main iteration
    [llh_,lambda_,eta_,delta_] = ...
        gh_llh(sample,mu,gamma,Sigma,lambda,chi,psi);
    llh(iter) = mean(llh_);
    
    if iter>2
        if abs(llh(iter)-llh(iter-1))<1e-5
            break;
        end
    end
    
    % a, b, c: n times 1
    [a,b,~] = gh_estep(lambda_,eta_,delta_);
    
    % sufficient statistics
    S1 = mean(a);
    S2 = mean(b);
    S5 = sample'*a/n;
    S6 = repmat(sqrt(a),1,dim).*sample;
    S6 = S6'*S6/n;
    
    % M step
    mu = (S4-S2*S5)/(1-S1*S2);
    gamma = (S5-S1*S4)/(1-S1*S2);
    Sigma = S6-S5*mu'-mu*S5'+S1*(mu*mu')-S2*(gamma*gamma');
    Sigma = setdet1(Sigma);
    
    % a, b, c: n times 1
    [a,b,c] = gh_estep(lambda_,eta_,delta_);
    
    % sufficient statistics
    S1 = mean(a);
    S2 = mean(b);
    S3 = mean(c);
    [lambda,chi,psi] = ...
        gig_mle(S1,S2,S3,[lambda;chi;psi],Aeq,beq,lb,ub,options);    
end

llh(iter+1:end) = llh(iter);

end
