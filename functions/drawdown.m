function [dd,loc] = drawdown(x)
n = length(x);
X = zeros(n,n);

for i = 1:n
    X(i:n,i) = x(1:(n-i+1));
end
X = tril(cumsum(X,2));
[dd,ind] = min(X);
loc = zeros(2,n);
loc(2,:) = ind;
loc(1,:) = ind-(0:n-1);
loc(:,dd>=0) = 0;
end

