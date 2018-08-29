%%
clear 
close all
clc
colors
% 
% %% Simulation parameters
% %% ------------------------------------------
% dt = 0.05;
% T = 40;%8*160;
% t = 0:dt:T;
% Steps = length(t);
% %% ------------------------------------------
% 
% %% Model parameters
% %% ------------------------------------------
% F = 8;
% n = 400;
% %% ------------------------------------------
% 
% %% Observations
% %% ------------------------------------------
% var_y = 1;
% skip = 2;
% H = getH(skip,n);
% R = var_y*eye(size(H,1));
% %% ------------------------------------------
% 
% %% Initial conditions/ spin up
% %% ------------------------------------------
% load LongSim400.mat
% %% ------------------------------------------
% 
% 
% % Generate data
% % ------------------------------------------
% Gap = 2;
% yo = yC(:,end);
% yAll = model(yo,dt,Steps,F);
% [z,tObs] = getObs(H,R,t,yAll,Gap,Steps);
% % %% ------------------------------------------
% 
% 
% %% EnKF 
% %% ------------------------------------------
% Ne = 100;
% locRad = 4;
% infl = 1.05;
% 
% Xo = yC(:,randi([1 length(yC)],1,Ne));
% Cloc = getCov(n,locRad); %GetLocMatrix2(n,locRad);
% 
% [x,traceP,X] = myEnKF(infl,Ne,Xo,Cloc,z,R,H,F,Gap,Steps,dt);
% MSE = mean((x(:,Gap+1:Gap:end) - yAll(:,Gap+1:Gap:end)).^2)';
% traceP = traceP/n;
% 
% % get rid of spin up
% SpinUp = 20;
%  
% fprintf('EnKF average MSE: %g\n',mean(MSE(SpinUp:end)))
% fprintf('EnKF average spread: %g\n',mean(traceP(SpinUp:end)))
% fprintf('EnKF average normalized MSE: %g\n',mean(MSE(SpinUp:end)./traceP(SpinUp:end)))
% disp(' ')
% 
% %
% figure(1)
% hold on, plot(t,yAll(1,:),'Color',Color(:,4),'LineWidth',2)
% hold on,plot(tObs,z(1,:),'.','Color',Color(:,4),'MarkerSize',20)
% hold on,plot(t,x(1,:),'-','Color',Color(:,2),'LineWidth',2)
% axis([t(SpinUp) t(end) -5 15])
% set(gcf,'Color','w'), box off, set(gca,'FontSize',20)
% 
% figure(2)
% hold on,plot(tObs(SpinUp:end),MSE(SpinUp:end),'Color',Color(:,1),'LineWidth',2)
% hold on,plot(tObs(SpinUp:end),traceP(SpinUp:end),':','Color',Color(:,1),'LineWidth',2)
% set(gcf,'Color','w'), box off, set(gca,'FontSize',20)
% %% ------------------------------------------
% 
% save EnKFRunVary1n40Gap2.mat

% load EnKFRunVary1n40Gap2.mat
% load EnKFRunVary1n400Gap2.mat
% load EnKFRunVary1n400Gap4.mat
% load EnKFRunVary1n400.mat
load EnKFRunVary1.mat
% load EnKFRunVary01.mat


%% Var PF
%% ------------------------------------------
yo = yAll(:,end);
yAll = model(yo,dt,Steps,F);
[z,tObs] = getObs(H,R,t,yAll,Gap,Steps);

Cloc = getCov(n,4);
nAssims = 100;
MSE = zeros(1,nAssims);
traceP = zeros(1,nAssims);

muo = x(:,end);
B = Cloc.*cov(X');
Lb = sqrtm(B);

%% Ne for EnKF in EnVar
Ne = 40;
X = X(:,1:Ne);
Cloc = getCov(n,4);
nObs = size(z,1);
for jj=1:nAssims 

    %% Var
    [mu,~,~,~,~,~,J] = myMinLS2(zeros(n,1),z(:,jj),Gap,dt,F,H,R,muo,Lb);
    Hess = 2*(J'*J);
    Cov = Hess\eye(n);
    Lcov = sqrt(1.0)*sqrtm(Cov);
    
    %% ensemble
    w = zeros(Ne,1);
    for kk=1:Ne
        xi = randn(n,1);
        tmp = mu + sqrt(1.0)*Lcov*xi;
        Fx = norm(funcF2(tmp,z(:,jj),Gap,dt,F,H,R,muo,Lb))^2;
        tmp = muo+Lb*tmp;
        Xo(:,kk) = tmp;
        tmp = model(tmp,dt,Gap+1,F);
        X(:,kk) = tmp(:,end);
        w(kk) = Fx - .5*norm(xi)^2;
    end
    w = normalizeweights(w);
%     XOLD = X;
%     X = resampling(w,X,Ne,n);   
    
    mu = muo+Lb*mu;
    Mmu = model(mu,dt,Gap+1,F);
    Mmu = Mmu(:,end);
        
    figure(20)
    hold off
    plot(1:2:n,z(:,jj),'.','Color',Color(:,4),'MarkerSize',20)
    hold on, plot(X,'Color',Color(:,2))
    hold on,plot(Mmu,'Color',Color(:,1),'LineWidth',2)
    hold on, plot(yAll(:,jj*Gap+1),'Color',Color(:,4))
    axis([1 n -10 15])
    drawnow
    
    
    %% local weights
    W = zeros(n,Ne);
    muP = muo+Lb*mu;
    CovP = Lb*Cov*Lb;
    Lp = sqrtm(CovP);
    for kk=1:n
        if mod(kk,2)~=0
            obsInd = mod(kk,2)*(kk+1)/2;
        else
            obsInd = kk/2;
        end
        locfx = QuadLoc(kk,n,1)';
        locfo = QuadLoc(obsInd,n/2,1)';
        
        wloc=zeros(Ne,1);
        for pp=1:Ne
             wloc(pp) = ...
               .5*norm(locfo.*(sqrt(R)\(z(:,jj)-H*X(:,pp))))^2 ...
                    + .5*norm( locfx.*(Lb\(muo -Xo(:,pp))) )^2 ...
                        -.5*norm( locfx.*(Lp\(muP -Xo(:,pp))) )^2;
        end
        wloc = normalizeweights(wloc);
        Rloc(kk)=mean(wloc.^2)/mean(wloc)^2;
        W(kk,:)=wloc';
    end
       
    U = zeros(n,Ne);
    for kk=1:Ne
        U(:,kk) = sqrt(W(:,kk)).*(X(:,kk)-mean(X,2));
    end
    
    %% update
    muo = Mmu;
    B =1.05*Cloc.*(U*U');
    Lb = sqrtm(B);

    %% Errors
    MSE(jj) = sum((muo-yAll(:,jj*Gap+1)).^2)/n;
    traceP(jj) = trace(B)/n;
    
    fprintf('MSE = %g\n',MSE(jj))
    fprintf('trac(P)/n = %g\n',traceP(jj) )
    fprintf('R =  %g\n',mean(w.^2)/mean(w)^2)
    fprintf('R loc = %g\n',mean(Rloc));
    disp(' ' )
    
end
figure
plot(MSE)
hold on, plot(traceP)
%%
disp(' ')
fprintf('Averages over %g assimilations\n',jj)
fprintf('MSE = %g\n',mean(MSE))
fprintf('trace(P)/n = %g\n',mean(traceP))
