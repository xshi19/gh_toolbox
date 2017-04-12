function [A,c] = setdet1(A)
% Set det(A)=1
d = size(A,1);
if det(A)==0
    c = det(A*eps^(-25/d))^(1/d)*eps^(25/d);
elseif det(A)==Inf
    c = det(A*eps^(25/d))^(1/d)*eps^(-25/d);
else
    c = det(A)^(1/d);
end
A = A/c;
end