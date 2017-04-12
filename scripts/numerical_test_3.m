% Test factor analysis with DJ30
clear; close all;
load('DJ30_2005_2015_RETURN.mat')

ret = data.return(end-2000:end,:);
ret = ret./(ones(size(ret,1),1)*std(ret));
dim = size(ret,2);
maxiter = 100;
numS = [500,500];

[mu,gamma,Sigma,lambda,chi,psi,llh] = gh_em(ret(:,1:dim));
sample = gh_rnd(sum(numS),mu,gamma,Sigma,lambda,chi,psi);

[mu1,gamma1,Sigma1,lambda1,chi1,psi1,llh1] = gh_em(sample(1:numS(1),:));
mu2 = mu1; gamma2 = gamma1; Sigma2 = Sigma1; 
lambda2 = lambda1; chi2 = chi1; psi2 = psi1;
l0 = sum(gh_llh(sample(1:numS(1),:),mu1,gamma1,Sigma1,lambda1,chi1,psi1));

t1 = 0;
t2 = 0;

for i = 1:numS(2)
    l1(i) = gh_llh(sample(i,:),mu1,gamma1,Sigma1,lambda1,chi1,psi1)-log(2*pi)*(dim/2);
    l2(i) = gh_llh(sample(i,:),mu2,gamma2,Sigma2,lambda2,chi2,psi2)-log(2*pi)*(dim/2);
    
    tic
    [mu1,gamma1,Sigma1,lambda1,chi1,psi1] = ...
        gh_em(sample(1:(numS(1)+i),:),mu1,gamma1,Sigma1,lambda1,chi1,psi1);
    t1 = t1+toc;
    
    tic
    [mu2,gamma2,Sigma2,lambda2,chi2,psi2] = ...
        gh_olem(sample(numS(1)+i,:),numS(1)+i,...
        mu2,gamma2,Sigma2,lambda2,chi2,psi2);
    t2 = t2+toc;
    
    lmle(i) = sum(gh_llh(sample(1:(numS(1)+i),:),...
        mu1,gamma1,Sigma1,lambda1,chi1,psi1));
    
    hd1(i) = gh_hdist(mu,gamma,Sigma,lambda,chi,psi,...
        mu1,gamma1,Sigma1,lambda1,chi1,psi1);
    hd2(i) = gh_hdist(mu,gamma,Sigma,lambda,chi,psi,...
        mu2,gamma2,Sigma2,lambda2,chi2,psi2);
end

figure
hold
plot(hd1)
plot(hd2,'r');
xlabel('Steps')
ylabel('Hellinger distance')
legend('EM', 'On-line EM')

cl1 = -cumsum(l1);
cl2 = -cumsum(l2);

figure
hold
plot(cl1)
plot(cl2,'r')
xlabel('Steps')
ylabel('Cumulative loss')
legend('EM', 'On-line EM')

results = [(50:50:500)',hd1(50:50:500)',...
    hd2(50:50:500)',cl1(50:50:500)',cl2(50:50:500)'];


