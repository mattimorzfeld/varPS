function [xAll,traceP,X]=myEnVar(infl,Cloc,Ne,X,muo,Lb,z,R,H,F,Gap,Steps,dt)
nAssims = size(z,2);
n = length(muo);
xAll = zeros(n,Steps-1);
traceP = zeros(nAssims,1);
for jj=1:nAssims 
     fprintf('Assim %g / %g\n',jj,nAssims);
    
    %% EnKF
    [~,~,X] = myEnKF(infl,Ne,X,Cloc,z(:,jj),R,H,F,Gap,Gap,dt);
    
    %% Var
    [mu,~,~,~,~,~,~] = myMinLS2(zeros(n,1),z(:,jj),Gap,dt,F,H,R,muo,Lb);    
    mu = muo+Lb*mu;
    Mmu = model(mu,dt,Gap+1,F);
    xAll(:,(jj-1)*Gap+1:jj*Gap+1)= Mmu;
    Mmu = Mmu(:,end);
     
    %% re-center ensemble
    X = X-mean(X,2)*ones(1,Ne);
    X = X+Mmu*ones(1,Ne);
        
    %% update
    muo = Mmu;
    P = infl*Cloc.*cov(X');
    Lb = real(sqrtm(P)');

    %% Save
    xAll(:,jj*Gap+1) = muo;
    traceP(jj) = trace(P)/n;
end
