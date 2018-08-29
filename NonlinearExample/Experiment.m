function [MSEw,Pw,MSEnw,Pnw,R,Skewness] = Experiment(Xt,r,Ne,e)

MSEw = zeros(length(Xt),1);
Pw = zeros(length(Xt),1);
MSEnw = zeros(length(Xt),1);
Pnw = zeros(length(Xt),1);
R = zeros(length(Xt),1);
Skewness = zeros(length(Xt),1);
for oo=1:length(Xt)
    xt = Xt(oo);
    y = M(xt,e)+1*sqrt(r)*randn;
    
    xr = -4:.001:4;
    Fr = F(xr,y,r,e);
    [a,b]=min(Fr);
    mu = xr(b);
    Hess = 1+(1/r)*(Mpp(mu,e)*(M(mu,e)-y) + Mp(mu,e)^2);
    C = 1/Hess;
    
    %% sampling
    Xs = zeros(Ne,1);
    w = zeros(Ne,1);
    for kk=1:Ne
        x = mu+sqrt(C)*randn;
        Xs(kk) = x;
        w(kk)  = F(x,y,r,e) - Fo(x,mu,C);
    end
    w = normalizeWeights(w);
    R(oo) = mean(w.^2)/mean(w)^2;
    Xrs = resampling(w,Xs',Ne,1);
    Xmw = Xs'*w;
    Pw(oo) = cov(Xrs);
    MSEw(oo) = (xt-Xmw)^2;
        
    Xmnw = mu;
    Pnw(oo) = cov(Xs);
    MSEnw(oo) = (xt-Xmnw)^2;

    Skewness(oo) = skewness(Xrs);

end
Pw = mean(Pw);
MSEw = mean(MSEw);
Pnw = mean(Pnw);
MSEnw = mean(MSEnw);
R = mean(R);
Skewness = mean(abs(Skewness));