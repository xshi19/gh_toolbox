clear; close all;
% load('DJ30_2005_2015_RETURN.mat')
load('sp500yahoo.mat')

d1 = datenum('2010-01-01');
d2 = datenum('2014-01-01');
d3 = datenum('2016-01-01');
r = 200;
winsize = 1000;

loc = logical((data.date<d3).*(data.date>=d1));
dates = data.date(loc);
sample = data.return(loc,:);
sample = sample(:,sum(isnan(sample))<=5);
sample(isnan(sample)) = 0;

[~,rep] = unique(sample(1,:));
sample = sample(:,sort(rep));
dim = size(sample,2);

tstart = find(data.date(loc)<d2,1,'last');

sample = sample./repmat(std(sample(1:tstart,:)),size(sample,1),1);
[mu1,gamma1,Sigma1,lambda1,chi1,psi1,llh] = gh_em(sample(1:tstart,:));
mu2 = mu1; gamma2 = gamma1; Sigma2 = Sigma1;
lambda2 = lambda1; chi2 = chi1; psi2 = psi1;
mu3 = mu1; gamma3 = gamma1; Sigma3 = Sigma1;
lambda3 = lambda1; chi3 = chi1; psi3 = psi1;

[V,L] = eigs(Sigma1,200);
F1 = V*sqrt(L);
D1 = diag(Sigma1-F1*F1');

[mu4,gamma4,Sigma4,lambda4,chi4,psi4,llh,F4,D4] = ...
    gh_fac(sample(1:tstart,:),mu1,gamma1,F1,D1,lambda1,chi1,psi1,50,'GH',0,1e-3);
mu5 = mu4; gamma5 = gamma4; Sigma5 = Sigma4; F5 = F4; D5 = D4;
lambda5 = lambda4; chi5 = chi4; psi5 = psi4;

loss = zeros(size(sample,1)-tstart,5);
time = zeros(size(sample,1)-tstart,5);
for t = 1:size(sample,1)-tstart
    datestr(dates(tstart+t))
    loss(t,1) = -gh_llh(sample(tstart+t,:),mu1,gamma1,Sigma1,lambda1,chi1,psi1);
    loss(t,2) = -gh_llh(sample(tstart+t,:),mu2,gamma2,Sigma2,lambda2,chi2,psi2);
    loss(t,3) = -gh_llh(sample(tstart+t,:),mu3,gamma3,Sigma3,lambda3,chi3,psi3);
    loss(t,4) = -gh_llh(sample(tstart+t,:),mu4,gamma4,Sigma4,lambda4,chi4,psi4);
    loss(t,5) = -gh_llh(sample(tstart+t,:),mu5,gamma5,Sigma5,lambda5,chi5,psi5);
    
    pastret = sample((tstart+t-winsize+1):(tstart+t),:);
    
    tic
    [mu1,gamma1,Sigma1,lambda1,chi1,psi1] = ...
        gh_em(pastret,mu1,gamma1,Sigma1,lambda1,chi1,psi1);
    time(t,1) = toc;
    
    tic
    [mu2,gamma2,Sigma2,lambda2,chi2,psi2] = ...
        gh_olem(sample(tstart+t,:),winsize+t,...
        mu2,gamma2,Sigma2,lambda2,chi2,psi2);
    time(t,2) = toc;
    
    tic
    [mu3,gamma3,Sigma3,lambda3,chi3,psi3] = ...
        gh_em(sample(1:tstart+t,:),mu3,gamma3,Sigma3,lambda3,chi3,psi3);
    time(t,3) = toc;
    
    tic
    [mu4,gamma4,Sigma4,lambda4,chi4,psi4,~,F4,D4] = ...
        gh_fac(pastret,mu4,gamma4,F4,D4,lambda4,chi4,psi4,50,'GH',0,1e-3);
    time(t,4) = toc;
    
    tic
    [mu5,gamma5,Sigma5,lambda5,chi5,psi5,~,F5,D5] = ...
        gh_fac(sample(1:tstart+t,:),mu5,gamma5,F5,D5,lambda5,chi5,psi5,50,'GH',0,1e-3);
    time(t,5) = toc;
    
end

cl = cumsum(loss);
dt = datetime(dates(tstart+1:end),'ConvertFrom','datenum');
y = unique(dt.Year);
m = unique(dt.Month);

result = [];
for i = 1:length(y);
    for j = 1:length(m)
        k = find((dt.Year==y(i)).*(dt.Month==m(j)),1,'last');
        if ~isempty(k)
            result = [result; [y(i),m(j),cl(k,:)]];
        end
    end
end

figure
hold
plot(dates(tstart+1:end),cl);
datetick('x')
ylabel('Cumulative loss')
legend('EM with moving window','On-line EM','EM with increasing sample',...
    'FA with moving window','FA with increasing sample')

