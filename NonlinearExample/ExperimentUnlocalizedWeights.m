function [MSEw,Pw,R] = ExperimentUnlocalizedWeights(Xt,r,Ne,e)
n = length(Xt);
w = zeros(Ne,1);
Xs = zeros(n,Ne);
for oo=1:n
    xt = Xt(oo);
    y = M(xt,e)+0*sqrt(r)*randn;
    
    xr = -6:.001:6;
    Fr = F(xr,y,r,e);
    [~,b]=min(Fr);
    mu = xr(b);
    Hess = 1+(1/r)*(Mpp(mu,e)*(M(mu,e)-y) + Mp(mu,e)^2);
    C = 1/Hess;
    
    %% sampling
    for kk=1:Ne
        x = mu+sqrt(C)*randn;
        Xs(oo,kk) = x;
        w(kk)  = w(kk) +F(x,y,r,e) - Fo(x,mu,C);
    end
end

w = normalizeWeights(w);
R = mean(w.^2)/mean(w)^2;
Xrs = resampling(w,Xs,Ne,n);
Xmw = Xs*w;%mean(Xrs,2);
Pw = trace(cov(Xrs'))/n;
MSEw = sum((Xt-Xmw).^2)/n;
        
