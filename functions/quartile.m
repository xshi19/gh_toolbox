function y = quartile(q,x)

xs = sort(x);
p = (1:length(x))'/length(x);
y = interp1(p,xs,q);

end