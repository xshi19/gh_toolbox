function [w,fval] = nm_cvaropt(w0,m,l,alpha,gamma,Sigma,y,z,c,A,b,Aeq,beq,lb,ub)


options = optimoptions('fmincon','GradObj','on',...
        'Display','off','Algorithm','interior-point');
[w,fval] = fmincon(@(w)cvarobj(w),w0,A,b,Aeq,beq,lb,ub,[],options);
if fval>=cvarobj(w0)
    disp('No change')
    w = w0;
end

function [f,g] = cvarobj(w)
[v,d] = nm_portcvar(w,alpha,gamma,Sigma,y,z);
f = -w'*m+l*v+c*norm(w-w0,1);
g = -m+l*d+c*sign(w-w0);
end

end