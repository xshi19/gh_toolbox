function [lambda,chi,psi,fval] = gig_mle(S1,S2,S3,x0,Aeq,beq,lb,ub,options)
% Maximum likelihood estimator for GIG
if nargin == 3
    Aeq = [];
    beq = [];
    lb = [-20;0;0];
    ub = [20;Inf;Inf];
    x0 = [0;1;1];
    options = optimoptions('fmincon','GradObj','on',...
        'Display','iter','Algorithm','interior-point');
end
fun = @(x)gig_llh(x(1),x(2),x(3),S1,S2,S3);
[xopt,fval] = fmincon(fun,x0,[],[],Aeq,beq,lb,ub,[],options);

if fval>gig_llh(x0(1),x0(2),x0(3),S1,S2,S3)
    disp('optimization incorrect')
    lambda = x0(1);
    chi = x0(2);
    psi = x0(3);
else
    lambda = xopt(1);
    chi = xopt(2);
    psi = xopt(3);
end

end