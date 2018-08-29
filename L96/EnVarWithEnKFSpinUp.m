clear 
close all
clc
colors

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
% n = 40;
% %% ------------------------------------------
% 
% %% Observations
% %% ------------------------------------------
% var_y = .1;
% skip = 2;
% H = getH(skip,n);
% R = var_y*eye(size(H,1));
% %% ------------------------------------------
% 
% %% Initial conditions/ spin up
% %% ------------------------------------------
% load LongSim.mat
% %% ------------------------------------------
% 
% 
% % Generate data
% % ------------------------------------------
% Gap = 1;
% yo = yC(:,end);
% yAll = model(yo,dt,Steps,F);
% [z,tObs] = getObs(H,R,t,yAll,Gap,Steps);
% % %% ------------------------------------------
% 
% 
% %% EnKF 
% %% ------------------------------------------
% Ne = 1000;
% locRad = 12;
% infl = 1.00;
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

% save EnKFRunVary01.mat


% load EnKFRunVary1n400Gap4.mat
% load EnKFRunVary1.mat
% load EnKFRunVary01.mat


% load EnKFRunVary1n40Gap2.mat
% load EnKFRunVary1n400Gap2.mat
load EnKFRunVary1n400Gap4.mat
% load EnKFRunVary1n400.mat
load EnKFRunVary1.mat
% load EnKFRunVary01.mat


%% Localized PF
%% ------------------------------------------
yo = yAll(:,end);
yAll = model(yo,dt,Steps,F);
[z,tObs] = getObs(H,R,t,yAll,Gap,Steps);

% plot(X,'b')
% hold on, plot(yo,'r')

Cloc = getCov(n,6);%GetLocMatrix2(n,40);
nAssims = 100;
MSE = zeros(1,nAssims);
traceP = zeros(1,nAssims);

muo = x(:,end);%yAll(:,1);%mean(X,2);
Lb = chol(Cloc.*cov(X'));

%% Ne for EnKF in EnVar
Ne = 20;
X = X(:,1:Ne);
Cloc = getCov(n,4);
for jj=1:nAssims 
    figure(20)
    hold off
    plot(1:2:n,z(:,jj),'.','Color',Color(:,4),'MarkerSize',20)
   
    %% EnKF
    [~,~,X] = myEnKF(1.05,Ne,X,Cloc,z(:,jj),R,H,F,Gap,Gap,dt);
    
    %% Var
    [mu,~,~,~,~,~,~] = myMinLS2(zeros(n,1),z(:,jj),Gap,dt,F,H,R,muo,Lb);    
    mu = muo+Lb*mu;
    Mmu = model(mu,dt,Gap+1,F);
    Mmu = Mmu(:,end);
    
     
    %% re-center ensemble
    X = X-mean(X,2)*ones(1,Ne);
    X = X+Mmu*ones(1,Ne);
    
    
    %% plot
    hold on, plot(X,'Color',Color(:,2))
    hold on,plot(Mmu,'Color',Color(:,1),'LineWidth',2)
    hold on, plot(yAll(:,jj*Gap+1),'Color',Color(:,4))
    axis([1 n -10 15])
    drawnow
    
    %% update
    muo = Mmu;
    P = Cloc.*cov(X');
    Lb = chol(P)';

    
    %% Errors
    MSE(jj) = sum((muo-yAll(:,jj*Gap+1)).^2)/n;
    traceP(jj) = trace(P)/n;
    
    fprintf('Step %g / %g\n',jj,nAssims)
    fprintf('MSE = %g\n',MSE(jj))
    fprintf('trac(P)/n = %g\n',traceP(jj) )
    disp(' ' )
end
figure
plot(MSE)
hold on, plot(traceP)
%%
disp(' ')
fprintf('Averages over %g assimilations\n',nAssims)
fprintf('MSE = %g\n',mean(MSE))
fprintf('trace(P)/n = %g\n',mean(traceP))
