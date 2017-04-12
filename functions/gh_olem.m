function [mu,gamma,Sigma,lambda,chi,psi] = ...
    gh_olem(sample,tau,mu,gamma,Sigma,lambda,chi,psi,distr,fixt)
% On-line (recursive) EM algorithm for GH distribution
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
if nargin<10
    fixt = 0;
end
if nargin<9
    distr = 'GH';
end

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

S4_ = mean(sample,1)';

[S1,S2,S3,S4,S5,S6] = param2sufstat(mu,gamma,Sigma,lambda,chi,psi);

[~,lambda_,eta_,delta_] = ...
        gh_llh(sample,mu,gamma,Sigma,lambda,chi,psi);

a = besselkratio(lambda_-1,lambda_,eta_);
b = besselkratio(lambda_+1,lambda_,eta_);
c = (besselkratio(lambda_+1e-10,lambda_,eta_)-...
    besselkratio(lambda_-1e-10,lambda_,eta_))/2e-10;
a = a./delta_;
b = b.*delta_;
c = c + log(delta_);

if c<Inf && c>-Inf
    S1_ = mean(a);
    S2_ = mean(b);
    S3_ = mean(c);
    S5_ = sample'*a/n;
    S6_ = repmat(sqrt(a),1,dim).*sample;
    S6_ = S6_'*S6_/n;
    
    S1 = S1 + (S1_-S1)/tau;
    S2 = S2 + (S2_-S2)/tau;
    S3 = S3 + (S3_-S3)/tau;
    S4 = S4 + (S4_-S4)/tau;
    S5 = S5 + (S5_-S5)/tau;
    S6 = S6 + (S6_-S6)/tau;
    
    mu = (S4-S2*S5)/(1-S1*S2);
    gamma = (S5-S1*S4)/(1-S1*S2);
    Sigma = S6-S5*mu'-mu*S5'+S1*(mu*mu')-S2*(gamma*gamma');
    
    if fixt==0
        try
            [lambda,chi,psi] = ...
                gig_mle(S1,S2,S3,[lambda;chi;psi],Aeq,beq,lb,ub,options);
        catch
            disp('Error in optimization.')
        end
    end
    [mu,gamma,Sigma,lambda,chi,psi] = ...
        gh_reg(mu,gamma,Sigma,lambda,chi,psi,4);
end
end



    
