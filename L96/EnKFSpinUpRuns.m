clear
close all
clc
colors

%% Simulation parameters
%% ------------------------------------------
dt = 0.05;
T = 100;
t = 0:dt:T;
Steps = length(t);
%% ------------------------------------------

%% Model parameters
%% ------------------------------------------
F = 8;
n = 40;
%% ------------------------------------------

%% Observations
%% ------------------------------------------
var_y = .1;
skip = 2;
H = getH(skip,n);
R = var_y*eye(size(H,1));
%% ------------------------------------------

%% Initial conditions/ spin up
%% ------------------------------------------
load(strcat('LongSim_n',num2str(n),'.mat'))
%% ------------------------------------------


% Generate data
% ------------------------------------------
Gap = 4;
yo = yC(:,end);
yAll = model(yo,dt,Steps,F);
[z,tObs] = getObs(H,R,t,yAll,Gap,Steps);
% %% ------------------------------------------


%% EnKF 
%% ------------------------------------------
Ne = 40;
locRad = 10.5;
infl = 1.12;

Xo = yC(:,randi([1 length(yC)],1,Ne));
Cloc = getCov(n,locRad); %GetLocMatrix2(n,locRad);

[x,traceP,X] = myEnKF(infl,Ne,Xo,Cloc,z,R,H,F,Gap,Steps,dt);
MSE = mean((x(:,Gap+1:Gap:end) - yAll(:,Gap+1:Gap:end)).^2)';
traceP = traceP/n;

% get rid of spin up
SpinUp = 20;
 
fprintf('EnKF average MSE: %g\n',mean(MSE(SpinUp:end)))
fprintf('EnKF average spread: %g\n',mean(traceP(SpinUp:end)))
fprintf('EnKF average normalized MSE: %g\n',mean(MSE(SpinUp:end)./traceP(SpinUp:end)))
disp(' ')

%
figure(1)
hold on, plot(t,yAll(1,:),'Color',Color(:,4),'LineWidth',2)
hold on,plot(tObs,z(1,:),'.','Color',Color(:,4),'MarkerSize',20)
hold on,plot(t,x(1,:),'-','Color',Color(:,2),'LineWidth',2)
axis([t(SpinUp) t(end) -5 15])
set(gcf,'Color','w'), box off, set(gca,'FontSize',20)

figure(2)
hold on,plot(tObs(SpinUp:end),MSE(SpinUp:end),'Color',Color(:,1),'LineWidth',2)
hold on,plot(tObs(SpinUp:end),traceP(SpinUp:end),':','Color',Color(:,1),'LineWidth',2)
set(gcf,'Color','w'), box off, set(gca,'FontSize',20)
%% ------------------------------------------

xo = yAll(:,end);
save(strcat('EnKFRunVary',num2str(var_y),'n',num2str(n),'Gap',num2str(Gap),'.mat'), 'X','xo')