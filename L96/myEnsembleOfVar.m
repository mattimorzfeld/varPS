function [xAll,traceP,X]=myEnsembleOfVar(infl,Cloc,Ne,X,muo,Lb,z,R,H,F,Gap,Steps,dt)
nAssims = size(z,2);
n = length(muo);
xAll = zeros(n,Steps-1);
traceP = zeros(nAssims,1);
nObs = size(z,1);
for jj=1:nAssims 
    tic
    fprintf('Assim %g / %g\n',jj,nAssims);
    

    %%
    Xv = zeros(n,Ne);
    for kk=2:Ne
        [muEn,~,~,~,~,~,~] = myMinLS2(zeros(n,1),z(:,jj)+sqrtm(R)*randn(nObs,1),...
                                            Gap,dt,F,H,R,X(:,kk), Lb);
        muEn = X(:,kk)+Lb*muEn;
        MmuEn = model(muEn,dt,Gap+1,F);
        MmuEn = MmuEn(:,end);
        Xv(:,kk) = MmuEn;
    end
        
    %% Var
    fprintf('starting optimization\n')
    [mu,~,~,~,~,~,~] = myMinLS2(zeros(n,1),z(:,jj),Gap,dt,F,H,R,muo, Lb);    
    mu = muo+Lb*mu;
    Mmu = model(mu,dt,Gap+1,F);
    xAll(:,(jj-1)*Gap+1:jj*Gap+1)= Mmu;
    Mmu = Mmu(:,end);
    Xv(:,1) = Mmu;
    
    %% EnKF
    [~,~,X] = myEnKF(infl,Ne,X,Cloc,z(:,jj),R,H,F,Gap,Gap,dt);
    
    %% re-center ensemble
    X = X-mean(X,2)*ones(1,Ne);
    X = X+Mmu*ones(1,Ne);


    
%     trace(cov(Xv')/n)
%     trace(cov(X')/n)
%     
%     figure(10)
%     hold off
%     plot(X,'b')
%     hold on, plot(Xv,'r')
%     hold on, plot(1:2:40,z(:,jj),'gx')
% 
%     if jj>1
        X = Xv;
%     end
%     pause
        
    %% update
    muo = Mmu;

    P = infl*Cloc.*cov(X');
    Lb = real(sqrtm(P)');

    %% Save
    xAll(:,jj*Gap+1) = muo;
    traceP(jj) = trace(P)/n;
    fprintf('Time for one assim: %g\n',toc)
end