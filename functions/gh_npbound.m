clear; close all;

d = 2;
n = 1e6;
rho = 5;

lambda = [-10, -0.5, -0.5, 0.5];
chi = [1, 1, 1, 1e-6];
psi = [1, 1e-6, 1, 1];

Sigma = eye(d);
Sigma(1,1) = rho;
Sigma(end,end) = 1/rho;


% generalized inverse Gaussian subordinators
T = ones(n, length(lambda));
pT = T;
for j = 1:length(lambda)
    % normalize T such that E[T]=1
    mT = gig_moment(1,lambda(j),chi(j),psi(j));
    chi(j) = chi(j)/mT;
    psi(j) = psi(j)*mT;
    
    T(:,j) = gig_rnd(n,lambda(j),chi(j),psi(j));
    mean(T(:,j))
    pT(:,j) = gig_pdf(T(:,j),lambda(j),chi(j),psi(j));
end

% plot the densities
figure
plot(T,pT,'.')
legend('thin tail','"inverse gamma"','"gamma"','inverse gaussian')

% normal rvs with different Sigma
X = randn(n,d);
X = [X; X*chol(Sigma)];


alpha = zeros(2*n,length(lambda));
beta = alpha;
for j = 1:length(lambda)
    % mixture the normal rvs
    Y = X.*repmat(sqrt(T(:,j)),2,d);
    
    % generalized hyperbolic pdfs
    p = gh_pdf(Y,zeros(d,1),zeros(d,1),eye(d),lambda(j),chi(j),psi(j));
    q = gh_pdf(Y,zeros(d,1),zeros(d,1),Sigma,lambda(j),chi(j),psi(j));
    
    [~,E] = sort(p./q);
    alpha(:,j) = cumsum(p(E)./(p(E)+q(E)))/n;
    beta(:,j)  = 1-cumsum(q(E)./(p(E)+q(E)))/n;
end

figure
plot(alpha, beta)
xlim([0,1])
ylim([0,1])
legend('thin tail','"inverse gamma"','"gamma"','inverse gaussian')








