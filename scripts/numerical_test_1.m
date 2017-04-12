% Test EM-algorithm with SP500 and FTSE Indices
clear; close all;
load('sp_ftse.mat')

ret = data.return(end-2000:end,:);
ret = ret./(ones(size(ret,1),1)*std(ret));
dim = size(ret,2);
maxiter = 100;
numS = 5000;

mu0 = mean(ret(:,1:dim))';
gamma0 = zeros(dim,1);
Sigma0 = cov(ret(:,1:dim));
lambda0 = (-10:10)';
chi0 = 1;
psi0 = 1;

llh = zeros(maxiter,length(lambda0));
llh1 = zeros(maxiter,length(lambda0));
ghhd1 = zeros(length(lambda0),1);
gighd1 = zeros(length(lambda0),1);
err1 = zeros(length(lambda0),6);

llh2 = zeros(maxiter,length(lambda0));
ghhd2 = zeros(length(lambda0),1);
gighd2 = zeros(length(lambda0),1);
err2 = zeros(length(lambda0),6);
t1 = 0; t2 = 0;

for i = 1:length(lambda0)
    i
    [mu,gamma,Sigma,lambda,chi,psi,llh(:,i)] = ...
        gh_em(ret(:,1:dim),mu0,gamma0,Sigma0,lambda0(i),chi0,psi0,maxiter,'GHl');
    sample = gh_rnd(numS,mu,gamma,Sigma,lambda,chi,psi);
    
    tic
    [mu1,gamma1,Sigma1,lambda1(i),chi1(i),psi1(i),llh1(:,i)] = gh_em(sample);
    t1 = t1+toc;
    
    tic
    [mu2,gamma2,Sigma2,lambda2(i),chi2(i),psi2(i),llh2(:,i)] = gh_mcecm(sample);
    t2 = t2+toc;
    
    % llh0(i) = mean(gh_llh(sample,mu,gamma,Sigma,lambda,chi,psi));
    % llh1(i) = mean(gh_llh(sample,mu1,gamma1,Sigma1,lambda1(i),chi1(i),psi1(i)));
    
    ghhd1(i) = gh_hdist(mu,gamma,Sigma,lambda,chi,psi,...
        mu1,gamma1,Sigma1,lambda1(i),chi1(i),psi1(i));
    gighd1(i) = gig_hdist(lambda,chi,psi,lambda1(i),chi1(i),psi1(i));
    ghhd2(i) = gh_hdist(mu,gamma,Sigma,lambda,chi,psi,...
        mu2,gamma2,Sigma2,lambda2(i),chi2(i),psi2(i));
    gighd2(i) = gig_hdist(lambda,chi,psi,lambda2(i),chi2(i),psi2(i));
    
    err1(i,1) = norm(mu1-mu)/norm(mu);
    err1(i,2) = norm(gamma1-gamma)/norm(gamma);
    err1(i,3) = norm(Sigma1-Sigma)/norm(Sigma);
    err1(i,4) = norm(lambda1(i)-lambda)/norm(lambda);
    err1(i,5) = norm(chi1(i)-chi)/norm(chi);
    err1(i,6) = norm(psi1(i)-psi)/norm(psi);
    
    err2(i,1) = norm(mu2-mu)/norm(mu);
    err2(i,2) = norm(gamma2-gamma)/norm(gamma);
    err2(i,3) = norm(Sigma2-Sigma)/norm(Sigma);
    err2(i,4) = norm(lambda2(i)-lambda)/norm(lambda);
    err2(i,5) = norm(chi2(i)-chi)/norm(chi);
    err2(i,6) = norm(psi2(i)-psi)/norm(psi);
end

w = [0;1];
X = sample*w;
[h,x] = hist(X,100);
h = h/sum(h);
p1 = gh_pdf(x',mu'*w,gamma'*w,w'*Sigma*w,lambda,chi,psi);
p2 = gh_pdf(x',mu1'*w,gamma1'*w,w'*Sigma1*w,lambda1(i),chi1(i),psi1(i));
p1 = p1/sum(p1);
p2 = p2/sum(p2);
figure
hold
bar(x,h)
plot(x,p1,'r-','LineWidth',2)
plot(x,p2,'g--','LineWidth',2)
legend('Sample histgram', 'True distribution', 'Fitted distribution')