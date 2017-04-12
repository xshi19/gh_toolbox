% Test factor analysis with DJ30
clear; close all;
load('DJ30_2005_2015_RETURN.mat')

ret = data.return(end-2000:end,:);
ret = ret./(ones(size(ret,1),1)*std(ret));
dim = size(ret,2);
maxiter = 100;
numS = 1000;

mu0 = mean(ret(:,1:dim))';
gamma0 = zeros(dim,1);
Sigma0 = cov(ret(:,1:dim));
lambda0 = (-10:10)';
chi0 = 1;
psi0 = 1;

llh = zeros(maxiter,length(lambda0));
llh1 = zeros(maxiter,length(lambda0));
llh2 = zeros(maxiter,length(lambda0));
llh3 = zeros(1,length(lambda0));

ghhd = zeros(length(lambda0),3);
gighd = zeros(length(lambda0),2);
condn = zeros(length(lambda0),3);

for i = 1:length(lambda0)
    i
    [mu,gamma,Sigma,lambda,chi,psi,llh(:,i)] = ...
        gh_em(ret(:,1:dim),mu0,gamma0,Sigma0,lambda0(i),chi0,psi0,maxiter,'GHl');
    [U,S,~] = svd(Sigma);
    S = 0.8*S+0.2*mean(diag(S))*eye(dim);
    Sigma = U*S*U';
    Sigma = setdet1(Sigma);
    sample = gh_rnd(numS,mu,gamma,Sigma,lambda,chi,psi);
    
    [mu1,gamma1,Sigma1,lambda1(i),chi1(i),psi1(i),llh1(:,i)] = gh_em(sample);
    [V,L] = eigs(Sigma1,10);
    F1 = V*sqrt(L);
    D1 = diag(Sigma1-F1*F1');
    [mu2,gamma2,Sigma2,lambda2(i),chi2(i),psi2(i),llh2(:,i)] = ...
        gh_fac(sample,mu1,gamma1,F1,D1,lambda1(i),chi1(i),psi1(i));
    
    % llh0(i) = mean(gh_llh(sample,mu,gamma,Sigma,lambda,chi,psi));
    % llh1(i) = mean(gh_llh(sample,mu1,gamma1,Sigma1,lambda1(i),chi1(i),psi1(i)));
    llh3(i) = mean(gh_llh(sample,mu1,gamma1,F1*F1'+diag(D1),lambda1(i),chi1(i),psi1(i)));
    ghhd(i,:) = [gh_hdist(mu,gamma,Sigma,lambda,chi,psi,...
        mu1,gamma1,Sigma1,lambda1(i),chi1(i),psi1(i)),...
        gh_hdist(mu,gamma,Sigma,lambda,chi,psi,...
        mu2,gamma2,Sigma2,lambda2(i),chi2(i),psi2(i)),...
        gh_hdist(mu,gamma,Sigma,lambda,chi,psi,...
        mu1,gamma1,F1*F1'+diag(D1),lambda1(i),chi1(i),psi1(i))];
    
    gighd(i,:) = [gig_hdist(lambda,chi,psi,lambda1(i),chi1(i),psi1(i)),...
        gig_hdist(lambda,chi,psi,lambda2(i),chi2(i),psi2(i))];
    
    condn(i,:) = [cond(Sigma1),cond(Sigma2),cond(F1*F1'+diag(D1))];

end

result = [lambda0,ghhd,condn,llh1(end,:)',llh2(end,:)',llh3'];
% w = [1;zeros(dim-1,1)];
w = ones(dim,1)/dim;
X = sample*w;
[h,x] = hist(X,100);
h = h/sum(h);
p0 = gh_pdf(x',mu'*w,gamma'*w,w'*Sigma*w,lambda,chi,psi);
p1 = gh_pdf(x',mu1'*w,gamma1'*w,w'*Sigma1*w,lambda1(i),chi1(i),psi1(i));
p2 = gh_pdf(x',mu2'*w,gamma2'*w,w'*Sigma2*w,lambda2(i),chi2(i),psi2(i));
p3 = gh_pdf(x',mu1'*w,gamma1'*w,w'*(F1*F1'+diag(D1))*w,lambda1(i),chi1(i),psi1(i));
p0 = p0/sum(p0);
p1 = p1/sum(p1);
p2 = p2/sum(p2);
p3 = p3/sum(p3);
figure
hold
bar(x,h)
plot(x,p0,'r-','LineWidth',2)
plot(x,p1,'g--','LineWidth',2)
plot(x,p2,'y:','LineWidth',2)
plot(x,p3,'m-.','LineWidth',2)
legend('Sample histgram', 'True distribution', 'EM', 'FA', 'PCA')