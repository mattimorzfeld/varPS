%%
clear 
close all
clc
colors

T =  5;
n = 40;
skip = 2;
var_y = .1;
Gap = 1;

% NeAll = 40;
% ClocRadAll = 7;
% LocRadAll = .1;
% inflAll = 1.04;
% Save = 0;
% EXIT = 0;
% [xAll,yAll,MSE,traceP,z,RLocSave,RSave] = OptimizeVarPF(n,T,var_y,skip,Gap,NeAll,LocRadAll,ClocRadAll,inflAll,Save,EXIT);


% NeAll = 30;
% ClocRadAll = [3];
% inflAll = [1.02];
% Save = 1;
% EXIT = 0;
% [xAll,yAll,MSE,traceP,z,RSave] = OptimizeVarPFNoLoc(n,T,var_y,skip,Gap,NeAll,ClocRadAll,inflAll,Save,EXIT);


% NeAll = 40;
% LocRadAll = 3;
% inflAll = 1.04;
% Save = 0;
% EXIT = 0;
% [xAll,yAll,MSE,traceP,z] = OptimizeEnKF(n,T,var_y,skip,Gap,NeAll,LocRadAll,inflAll,Save,EXIT);

% NeAll = [10 40];
% LocRadAll = [3.4 3.6];
% inflAll = [1.04 1.06];
% Save = 1;
% EXIT = 0;
% [xAll,yAll,MSE,traceP,z] = OptimizeEnVar(n,T,var_y,skip,Gap,NeAll,LocRadAll,inflAll,Save,EXIT);


% NeAll = [10 40];
% LocRadAll = [4];
% inflAll = [0.4];
% Save = 1;
% EXIT = 0;
% [xAll,yAll,MSE,traceP,z] = OptimizePoterjoysPF(T,var_y,skip,Gap,NeAll,LocRadAll,inflAll,Save,EXIT);


NeAll = [10 40];
LocRadAll = [4];
inflAll = 0.4;
Save = 1;
EXIT = 0;
[xAll,yAll,MSE,traceP,z] = OptimizePoterjoysPF_KDDM(T,var_y,skip,Gap,NeAll,LocRadAll,inflAll,Save,EXIT);

%%
t = 0:.05:T;
tObs = t(Gap+1:Gap:end);
figure(1)
subplot(211)
hold on, plot(t,yAll(1,:),'Color',Color(:,4),'LineWidth',2)
hold on,plot(tObs,z(1,:),'.','Color',Color(:,4),'MarkerSize',20)
hold on,plot(t,xAll(1,:),'-','Color',Color(:,2),'LineWidth',2)
subplot(212)
hold on, plot(t,yAll(2,:),'Color',Color(:,4),'LineWidth',2)
hold on,plot(t,xAll(2,:),'-','Color',Color(:,2),'LineWidth',2)
set(gcf,'Color','w')
set(gca,'FontSize',20)
box off

figure
plot(tObs,MSE,'Color',Color(:,2),'LineWidth',2)
hold on, plot(tObs,traceP,'--','Color',Color(:,1),'LineWidth',2)
set(gcf,'Color','w')
set(gca,'FontSize',20)
box off

disp(' ')
fprintf('MSE = %g\n',mean(MSE))
fprintf('trace(P)/n = %g\n',mean(traceP))
set(gcf,'Color','w')
set(gca,'FontSize',20)
box off

