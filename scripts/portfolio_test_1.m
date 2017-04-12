clear; close all;
load('DJ30_2005_2015_RETURN.mat')

d1 = datenum('2010-01-01');
d2 = datenum('2014-01-01');
d3 = datenum('2015-01-01');
winsize = 500;
l = 0.1;
c = 0.01;
alpha = 0.01;

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
[mu,gamma,Sigma,lambda,chi,psi,llh] = gh_em(sample(1:tstart,:));

lb = zeros(dim,1);
ub = ones(dim,1);
Aeq = ones(1,dim);
beq = 1;

% T = 100;
T = size(sample,1)-tstart-1;
Z = normrnd(0,1,1e5,1);
w1 = ones(dim,T+1)/dim;
w2 = ones(dim,T+1)/dim;
er = zeros(T,3);
cv = zeros(T,3);
to = zeros(T,3);
fv = zeros(T,3);
t1 = 0;
t2 = 0;

    
for t = 1:T
    % pastret = sample((tstart+t-winsize+1):(tstart+t),:);
    t
    w01 = w1(:,t).*exp(sample(tstart+t,:)');
    w02 = w1(:,t).*exp(sample(tstart+t,:)');
    w01 = w01/sum(w01);
    w02 = w02/sum(w02);
    
    [mu,gamma,Sigma,lambda,chi,psi] = ...
        gh_olem(sample(tstart+t,:),winsize+t,...
        mu,gamma,Sigma,lambda,chi,psi);
    
    Y = gig_rnd(1e5,lambda,chi,psi);
    [v,d,H] = nm_portcvar(w01,alpha,gamma,Sigma,Y,Z);
    m = gh_meanvar(mu,gamma,Sigma,lambda,chi,psi);
    m = m+l*mu;
    
    tic
    w1(:,t+1) = portprog(w01,m-l*d,l*H,c,[],[],Aeq,lb,ub);
    t1 = t1+toc
    tic
    w2(:,t+1) = nm_cvaropt(w02,m,l,alpha,...
        gamma,Sigma,Y,Z,c,[],[],Aeq,1,lb,ub);
    t2 = t2+toc
    er(t,1) = m'*w1(:,t+1);
    cv(t,1) = nm_portcvar(w1(:,t+1),alpha,gamma,Sigma,Y,Z);
    to(t,1) = norm(w1(:,t+1)-w01,1);
    fv(t,1) = er(t,1)-l*cv(t,1)-c*to(t,1);
    
    er(t,2) = m'*w2(:,t+1);
    cv(t,2) = nm_portcvar(w2(:,t+1),alpha,gamma,Sigma,Y,Z);
    to(t,2) = norm(w2(:,t+1)-w02,1);
    fv(t,2) = er(t,2)-l*cv(t,2)-c*to(t,2);
    
    er(t,3) = m'*w1(:,1);
    cv(t,3) = nm_portcvar(w1(:,1),alpha,gamma,Sigma,Y,Z);
    fv(t,3) = er(t,3)-l*cv(t,3)-c*to(t,3);
end


portret = [sum(sample(tstart+(1:T+1),:).*w1(:,1:T+1)',2),...
    sum(sample(tstart+(1:T+1),:).*w2(:,1:T+1)',2),...
    mean(sample(tstart+(1:T+1),:),2)];

figure
subplot(2,2,1)
plot(dates(tstart+(1:T)),er)
datetick('x','mmm-yy')
ylabel('Expected return')
legend('Approximation','True','Equal weight')

subplot(2,2,2)
plot(dates(tstart+(1:T)),cv)
datetick('x','mmm-yy')
ylabel('Expected CVaR')
legend('Approximation','True','Equal weight')

subplot(2,2,3)
plot(dates(tstart+(1:T)),to)
datetick('x','mmm-yy')
ylabel('Turnover')
legend('Approximation','True','Equal weight')

subplot(2,2,4)
plot(dates(tstart+(1:T)),fv)
datetick('x','mmm-yy')
ylabel('Optimized function value')
legend('Approximation','True','Equal weight')

figure
plot(dates(tstart+(1:T+1)),cumsum(portret(:,1:3)))
datetick('x','mmm-yy')
legend('Approximation','True','Equal weight')

figure
h = stairs([w1(:,end),w2(:,end)]);
h(1).Marker = 'o';
h(2).Marker = '*';
legend('Approximation','True')
xlabel('Assets')
ylabel('Weights')



