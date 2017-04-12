function [w,fval] = portprog(w0,mu,Sigma,c,A,b,Aeq,lb,ub)
% Solve min -mu'*w + 1/2*w'*Sigma*w + c*norm(w-w0,1)
% s.t. A*w<=b; Aeq*w==Aeq*w0

d = length(mu);
mu_ = -mu+Sigma*w0;
mu_ = [mu_;-mu_]+c;
Sigma_ = [Sigma,-Sigma;-Sigma,Sigma];
% [U,S,~] = svd(Sigma_);
% S = diag(S);
% S(S<1e-20) = 1e-20;
% Sigma_ = U*diag(S)*U';

A = [A;eye(d);-eye(d)];
b = [b;ub;-lb];
A_ = [A,-A];
b_ = b-A*w0;
Aeq_ = [Aeq,-Aeq];
beq_ = zeros(size(Aeq_,1),1);

if any(b<A*w0)
    disp('w0 does not satisfy the constraints')
end

opts = optimoptions('quadprog','Algorithm',...
    'interior-point-convex', 'Display','off');
[x,v] = quadprog(Sigma_,mu_,A_,b_,Aeq_,beq_,zeros(2*d,1),[],zeros(2*d,1),opts);

% figure
% hold
% bar(x(1:d))
% bar(-x(d+1:2*d))

if v>0
    disp('No change')
    w = w0;
else
    w = x(1:d)-x(d+1:2*d)+w0;
end
fval = -mu'*w+1/2*w'*Sigma*w+c*norm(w-w0,1);

end



