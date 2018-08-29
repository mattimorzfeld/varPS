function [xAll,traceP,X,RSave,RLocSave]=myVarPF(infl,locRad,Cloc,Ne,muo,Lb,z,R,H,F,Gap,Steps,dt,skip)
nObs = size(z,1);
nAssims = size(z,2);
n = length(muo);
xAll = zeros(n,Steps-1);
RLocSave = zeros(nAssims,1);
RSave = zeros(nAssims,1);
traceP = zeros(nAssims,1);
for jj=1:nAssims 
    fprintf('Assim %g / %g\n',jj,nAssims);
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
    
    %% CAREFUL WORKS FOR skip = 2 ONLY
%     disp('Starting weight computation')
    muP = muo+Lb*mu;
    CovP = Lb*Cov*Lb;
    Lp = sqrtm(CovP);
    W = zeros(nObs,Ne);
    Rloc = zeros(nObs,1);
    for oo=1:nObs
        wloc=zeros(Ne,1);
        
        ind = (oo-1)*skip+1;
        locf = QuadLoc(ind,n,locRad)';
       for kk=1:Ne      
           wloc(kk) = ...
               .5*norm(sqrt(R(oo,oo))\(z(oo,jj)-H(oo,:)*X(:,kk)))^2 ...
                    +.5*norm( locf.*(Lb\(muo -Xo(:,kk))) )^2 ...
                        -.5*norm( locf.*(Lp\(muP -Xo(:,kk))) )^2;
       end
       wloc = normalizeweights(wloc);
       Rloc(oo)=mean(wloc.^2)/mean(wloc)^2;
       W(oo,:)=wloc';
    end
    
    Wx = zeros(n,Ne);
    for kk=1:nObs
        Wx((kk-1)*skip+1:kk*skip,:) = [W(kk,:);W(kk,:)];
    end

    mx  = zeros(n,1);
    for aa=1:Ne
        mx = mx+Wx(:,aa).*X(:,aa);
    end

    U = zeros(n,Ne);
    for kk=1:Ne
        U(:,kk) = sqrt(Wx(:,kk)).*(X(:,kk)-mx);
    end
%     disp('End weight computation')

    mu = muo+Lb*mu;
    Mmu = model(mu,dt,Gap+1,F);
    Mmu = Mmu(:,end);
        
    
    %% update
    muo = Mmu;
    B = infl*Cloc.*(U*U');
    Lb = sqrtm(B);

    
    %% Save
%     MSE(jj) = sum((muo-yAll(:,jj*Gap+1)).^2)/n;
    xAll(:,jj*Gap+1) = muo;
    traceP(jj) = trace(B)/n;
    RSave(jj) = mean(w.^2)/mean(w)^2;
    RLocSave(jj) = mean(Rloc);

    
%     fprintf('MSE = %g\n',MSE(jj))
%     fprintf('trac(P)/n = %g\n',traceP(jj) )
%     fprintf('R loc = %g\n',mean(Rloc));
%     fprintf('R =  %g\n',mean(w.^2)/mean(w)^2)
%     disp(' ' )

	save('varPF2000Result.mat')    
end
