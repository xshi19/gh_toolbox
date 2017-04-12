clear; close all;
load('DJ30_2005_2015_RETURN.mat')
% load('sp500yahoo.mat')
% 
% d1 = datenum('2005-01-01');
% d2 = datenum('2008-01-01');
% d3 = datenum('2011-01-01');

% d1 = datenum('2005-01-01');
% d2 = datenum('2008-01-01');
% d3 = datenum('2011-01-01');
d1 = datenum('2008-01-01');
d2 = datenum('2011-01-01');
d3 = datenum('2014-01-01');
winsize = 500;
alpha = 0.01;

loc = logical((data.date<d3).*(data.date>=d1));
dates = data.date(loc);
sample = data.return(loc,:);
sample = sample(:,sum(isnan(sample))<=5);
sample(isnan(sample)) = 0;

[~,rep] = unique(sample(1,:));
sample = sample(:,sort(rep));
dim = size(sample,2);
% dim = 100;
% sample = sample(:,101:dim+100);

tstart = find(data.date(loc)<d2,1,'last');

sample_ = sample./repmat(std(sample(1:tstart,:)),size(sample,1),1);

[~,~,~,lambda,chi,psi] = gh_em(mean(sample_(1:tstart,:),2));
mu0 = mean(sample_(1:tstart,:))';
gamma0 = zeros(dim,1);
Sigma0 = cov(sample_(1:tstart,:))/gig_moment(1,lambda,chi,psi);
[mu,gamma,Sigma,lambda,chi,psi] = gh_em(sample_(1:tstart,:),...
    mu0,gamma0,Sigma0,lambda,chi,psi,100,'GH',1);
mu = mu.*std(sample(1:tstart,:))';
gamma = gamma.*std(sample(1:tstart,:))';
Sigma = diag(std(sample(1:tstart,:)))*Sigma*diag(std(sample(1:tstart,:)));
[m,C] = gh_meanvar(mu,gamma,Sigma,lambda,chi,psi);

lb = zeros(dim,1);
ub = ones(dim,1);
Aeq = [ones(1,dim);m'];
beq = [1;mean(m)];

Y = gig_rnd(1e5,lambda,chi,psi);
Z = normrnd(0,1,1e5,1);
options = optimoptions('fmincon','GradObj','off',...
        'Display','iter','Algorithm','interior-point');
w0 = ones(dim,1)/dim;
w_enb = fmincon(@(w)-enbCov(w,C),w0,[],[],Aeq,beq,lb,ub,[],options);
w_cvar = fmincon(@(w)nm_portcvar(w,alpha,gamma,Sigma,Y,Z)-w'*mu,w0,[],[],Aeq,beq,lb,ub,[],options);
w_cvarenb = fmincon(@(w)-nm_cvarenb(w,alpha,mu,gamma,Sigma,Y,Z),w0,[],[],Aeq,beq,lb,ub,[],options);

[y,x] = hist(sample(1:tstart,:)*w0,200);
p = gh_pdf(x',mu'*w0,gamma'*w0,w0'*Sigma*w0,lambda,chi,psi);
y = y/sum(y);
p = p/sum(p);
figure
hold
stairs(x,y)
plot(x,p,'r')

results = [enbCov(w0,C),enbCov(w_enb,C),...
    enbCov(w_cvar,C),enbCov(w_cvarenb,C);...
    nm_cvarenb(w0,alpha,mu,gamma,Sigma,Y,Z),...
    nm_cvarenb(w_enb,alpha,mu,gamma,Sigma,Y,Z),...
    nm_cvarenb(w_cvar,alpha,mu,gamma,Sigma,Y,Z),...
    nm_cvarenb(w_cvarenb,alpha,mu,gamma,Sigma,Y,Z);...
    nm_portcvar(w0,alpha,gamma,Sigma,Y,Z)-w0'*mu,...
    nm_portcvar(w_enb,alpha,gamma,Sigma,Y,Z)-w_enb'*mu,...
    nm_portcvar(w_cvar,alpha,gamma,Sigma,Y,Z)-w_cvar'*mu,...
    nm_portcvar(w_cvarenb,alpha,gamma,Sigma,Y,Z)-w_cvarenb'*mu]'

% ret_in = sample(1:tstart,:)*[w0,w_enb,w_cvar,w_cvarenb];
% ret_out = sample(tstart+1:end,:)*[w0,w_enb,w_cvar,w_cvarenb];
ret_in = log(1+(exp(sample(1:tstart,:))-1)*[w0,w_enb,w_cvar,w_cvarenb]);
ret_out = log(1+(exp(sample(tstart+1:end,:))-1)*[w0,w_enb,w_cvar,w_cvarenb]);
dd_in = zeros(4,size(ret_in,1));
dd_out = zeros(4,size(ret_out,1));
for i = 1:4
    dd_in(i,:) = drawdown(ret_in(:,i));
    dd_out(i,:) = drawdown(ret_out(:,i));
end

% ret_in_mv = ret_in;
% ret_out_mv = ret_out;
% 
% for i = 1:4
%     ret_in_mv(:,i) = smooth(ret_in(:,i),5);
%     ret_out_mv(:,i) = smooth(ret_out(:,i),5);
% end

q = (0.01:0.001:0.1)';
q_in = quartile(q,dd_in');
q_out = quartile(q,dd_out');
q_in = [min(dd_in'); q_in];
q_out = [min(dd_out'); q_out];

% q = (0.01:0.001:0.1)';
% q_in = quartile(q,ret_in);
% q_out = quartile(q,ret_out);
% q_in = [min(ret_in); q_in];
% q_out = [min(ret_out); q_out];

figure
hold
plot([0;q],-q_in);
xlim([0,max(q)])
legend('Equally weighted','Best variance-based ENB','Smallest CVaR', 'Best CVaR-based ENB') 
set(gca,'xticklabel',cellstr(num2str(get(gca,'xtick')'*100,'%.1f%%')))
xlabel('Percentage')
ylabel('Drawdown')

figure
plot([0;q],-q_out);
xlim([0,max(q)])
legend('Equally weighted','Best variance-based ENB','Smallest CVaR', 'Best CVaR-based ENB') 
set(gca,'xticklabel',cellstr(num2str(get(gca,'xtick')'*100,'%.1f%%')))
xlabel('Percentage')
ylabel('Drawdown')

figure
hold
stairs([w0,w_enb,w_cvar,w_cvarenb])
xlim([1,dim])
legend('Equally weighted','Best variance-based ENB','Smallest CVaR', 'Best CVaR-based ENB') 
xlabel('Stock')
ylabel('Weight')

figure
plot(cumsum(ret_in))

figure
plot(cumsum(ret_out))

