function [xAll,traceP,X,RSave]=myVarPFNoLoc(infl,Cloc,Ne,muo,Lb,z,R,H,F,Gap,Steps,dt)
nAssims = size(z,2);
n = length(muo);
xAll = zeros(n,Steps-1);
RSave = zeros(nAssims,1);
traceP = zeros(nAssims,1);
for jj=1:nAssims 
%     fprintf('Assim %g / %g\n',jj,nAssims);
    %% Var
    [mu,~,~,~,~,~,J] = myMinLS2(zeros(n,1),z(:,jj),Gap,dt,F,H,R,muo,Lb);
    Hess = 2*(J'*J);
    Cov = Hess\eye(n);
    Lcov = sqrt(1.0)*sqrtm(Cov);
    
    %% ensemble
    w = zeros(Ne,1);
    Xo = zeros(n,Ne);
    X = zeros(n,Ne);
    RunningMean = zeros(n,Gap+1);
    for kk=1:Ne
        xi = randn(n,1);
        tmp = mu + sqrt(1.0)*Lcov*xi;
        Fx = norm(funcF2(tmp,z(:,jj),Gap,dt,F,H,R,muo,Lb))^2;
        tmp = muo+Lb*tmp;
        Xo(:,kk) = tmp;
        tmp = model(tmp,dt,Gap+1,F);
        RunningMean=RunningMean+tmp;
        X(:,kk) = tmp(:,end);
        w(kk) = Fx - .5*norm(xi)^2;
    end
    RunningMean = RunningMean/Ne;
    xAll(:,(jj-1)*Gap+1:jj*Gap+1)=RunningMean;
    w = normalizeweights(w);    

    mu = muo+Lb*mu;
    Mmu = model(mu,dt,Gap+1,F);
    Mmu = Mmu(:,end);
        
    
    %% update
    X = resampling(w,X,Ne,n);
    muo = Mmu;
    B = infl*Cloc.*cov(X');
    Lb = real(sqrtm(B));

    
    %% Save
    xAll(:,jj*Gap+1) = muo;
    traceP(jj) = trace(B)/n;
    RSave(jj) = mean(w.^2)/mean(w)^2;

end
