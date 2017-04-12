clear; close all;
% load('DJ30_2005_2015_RETURN.mat')
load('sp500yahoo.mat')

d1 = datenum('2010-01-01');
d2 = datenum('2014-01-01');
d3 = datenum('2016-01-01');

loc = logical((data.date<d3).*(data.date>=d1));
dates = data.date(loc);
sample = data.return(loc,:);
sample = sample(:,sum(isnan(sample))<=5);
sample(isnan(sample)) = 0;

[~,rep] = unique(sample(1,:));
sample = sample(:,sort(rep));
dim = size(sample,2);

sample_in = sample(data.date(loc)<d2,:);
sample_out = sample(data.date(loc)>=d2,:);

sample_out = sample_out./(ones(size(sample_out,1),1)*std(sample_in));
sample_in = sample_in./(ones(size(sample_in,1),1)*std(sample_in));

maxiter = 100;
mu0 = mean(sample_in)';
gamma0 = zeros(dim,1);
lambda0 = -0.5;
chi0 = 1;
psi0 = 1;
Sigma0 = cov(sample_in)/gig_moment(1,lambda0,chi0,psi0);

[mu1,gamma1,Sigma1,lambda1,chi1,psi1,llh1] = ...
    gh_em(sample_in,mu0,gamma0,Sigma0,lambda0,chi0,psi0,maxiter,'GH');

[mu2,gamma2,Sigma2,lambda2,chi2,psi2,llh2] = ...
    gh_em(sample_in,mu0,gamma0,Sigma0,lambda0,chi0,psi0,maxiter,'NIG');

[mu3,gamma3,Sigma3,lambda3,chi3,psi3,llh3] = ...
    gh_em(sample_in,mu0,gamma0,Sigma0,lambda0,chi0,psi0,maxiter,'VG');

figure
plot([llh1,llh2,llh3])

w = unifrnd(-1,1,100,dim);
w = w./repmat(sum(w,2),1,dim);
Z = randn(1e6,1);
Y1 = gig_rnd(1e6,lambda1,chi1,psi1);
Y2 = gig_rnd(1e6,lambda2,chi2,psi2);
Y3 = gig_rnd(1e6,lambda3,chi3,psi3);

for i = 1:100
    i
    X1 = w(i,:)*mu1+w(i,:)*gamma1*Y1+sqrt(w(i,:)*Sigma1*w(i,:)')*sqrt(Y1).*Z;
    X2 = w(i,:)*mu2+w(i,:)*gamma2*Y2+sqrt(w(i,:)*Sigma2*w(i,:)')*sqrt(Y2).*Z;
    X3 = w(i,:)*mu3+w(i,:)*gamma3*Y3+sqrt(w(i,:)*Sigma3*w(i,:)')*sqrt(Y3).*Z;
    X4 = mean(sample_in*w(i,:)')+std(sample_in*w(i,:)')*Z;
    
    [ksH_in(i,1),pval_in(i,1),ksstat_in(i,1)] = kstest2(sample_in*w(i,:)',X1);
    [ksH_in(i,2),pval_in(i,2),ksstat_in(i,2)] = kstest2(sample_in*w(i,:)',X2);
    [ksH_in(i,3),pval_in(i,3),ksstat_in(i,3)] = kstest2(sample_in*w(i,:)',X3);
    [ksH_in(i,4),pval_in(i,4),ksstat_in(i,4)] = kstest2(sample_in*w(i,:)',X4);
    
    [ksH_out(i,1),pval_out(i,1),ksstat_out(i,1)] = kstest2(sample_out*w(i,:)',X1);
    [ksH_out(i,2),pval_out(i,2),ksstat_out(i,2)] = kstest2(sample_out*w(i,:)',X2);
    [ksH_out(i,3),pval_out(i,3),ksstat_out(i,3)] = kstest2(sample_out*w(i,:)',X3);
    [ksH_out(i,4),pval_out(i,4),ksstat_out(i,4)] = kstest2(sample_out*w(i,:)',X4);
    
    [adH_in(i,1),adstat_in(i,1)] = adtest2(sample_in*w(i,:)',X1);
    [adH_in(i,2),adstat_in(i,2)] = adtest2(sample_in*w(i,:)',X2);
    [adH_in(i,3),adstat_in(i,3)] = adtest2(sample_in*w(i,:)',X3);
    [adH_in(i,4),adstat_in(i,4)] = adtest2(sample_in*w(i,:)',X4);
    
    [adH_out(i,1),adstat_out(i,1)] = adtest2(sample_out*w(i,:)',X1);
    [adH_out(i,2),adstat_out(i,2)] = adtest2(sample_out*w(i,:)',X2);
    [adH_out(i,3),adstat_out(i,3)] = adtest2(sample_out*w(i,:)',X3);
    [adH_out(i,4),adstat_out(i,4)] = adtest2(sample_out*w(i,:)',X4);
end

result1 = [sum(ksH_in)',mean(pval_in)',mean(ksstat_in)',...
    sum(ksH_out)',mean(pval_out)',mean(ksstat_out)']

result2 = [sum(adH_in)',mean(adstat_in)',...
    sum(adH_out)',mean(adstat_out)']

param = [lambda1,chi1,psi1;lambda2,chi2,psi2;lambda3,chi3,psi3;]

figure
[y,x] = hist(sample_in*w(2,:)',100);
p1 = gh_pdf(x',w(2,:)*mu1,w(2,:)*gamma1,w(2,:)*Sigma1*w(2,:)',lambda1,chi1,psi1);
p2 = gh_pdf(x',w(2,:)*mu2,w(1,:)*gamma2,w(2,:)*Sigma2*w(2,:)',lambda2,chi2,psi2);
p3 = gh_pdf(x',w(2,:)*mu3,w(1,:)*gamma3,w(2,:)*Sigma3*w(2,:)',lambda3,chi3,psi3);
p4 = normpdf(x',mean(sample_in*w(2,:)'),std(sample_in*w(2,:)'));
y = y/sum(y);
p1 = p1/sum(p1);
p2 = p2/sum(p2);
p3 = p3/sum(p3);
p4 = p4/sum(p4);
p = [p1,p2,p3,p4];

figure
hold
bar(x,y)
plot(x,p,'LineWidth',1.5)
legend('sample','GH','NIG','VG','Gaussian')





