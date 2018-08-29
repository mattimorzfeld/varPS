function [Xo,yo] = SpinUp(NeSpinUp,TSpinUp,dt,F,stoch,Lq,H,R,Gap,n,Q,Plot)
colors
%% Spin-up: Use EnKF to generate initial ensemble
%% ------------------------------------------
load LongSim.mat
SpinUpSteps = TSpinUp/dt+1;
tSpinUp = 0:dt:TSpinUp;
yAll = model(yC(:,end),dt,SpinUpSteps,F,stoch,Lq);
[z,tObs] = getObs(H,R,tSpinUp,yAll,Gap,SpinUpSteps); % create synthetic data

Xo = zeros(n,NeSpinUp);
for kk=1:NeSpinUp
    jj = floor(1+length(yC)*rand);
    Xo(:,kk) = yC(:,jj);
end

Cloc = GetLocMatrix(n,.8,2);
LocRad = 0;

[xEnKF,traceP,X,minR,maxR,meanR] = myEnKFCompWeights(LocRad,NeSpinUp,Xo,Cloc,z,R,H,F,Gap,SpinUpSteps,dt,Q,stoch);
if Plot == 1
    figure(9)
    plot(tSpinUp,yAll(1,:),'Color',Color(:,2),'LineWidth',2)
    hold on,plot(tObs,z(1,:),'.','Color',Color(:,4),'MarkerSize',20)
    hold on,plot(tSpinUp,xEnKF(1,:),'-','Color',Color(:,4),'MarkerSize',20)
    set(gcf,'Color','w'), box off, set(gca,'FontSize',20)
    
    MSE = mean((xEnKF(:,Gap+1:Gap:end) - yAll(:,Gap+1:Gap:end)).^2);
    traceP = traceP'/n;
    figure
    plot(MSE./traceP)
    mean(MSE(end/2:end)./traceP(end/2:end))
end
yo = yAll(:,end);

%% ------------------------------------------