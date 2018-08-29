%%
clear
close all
clc
colors

n = 100;

%% The system
A = speye(n);
Q = speye(n);
R = 0.1*speye(n);
H = speye(n);
k = size(H,1); % number of obs

%% KF
nTries  = 1000; % compute MSE and traceP nTries times to average

MSESave= zeros(nTries,1);
tracePSave= zeros(nTries,1);
MSEdTraceP = zeros(nTries,1);
for kk=1:nTries
    % "Truth"
    xo = randn(n,1);
    xt = A*xo + Q*randn(n,1);
    y =  H*xt+sqrt(R)*randn(n,1);
    % filter
    [Xam,traceP] = myKF(y,A,Q,H,R);
    MSE = sum((Xam - xt).^2)/n;
    
    MSESave(kk) = MSE;
    tracePSave(kk) = traceP;    
    MSEdTraceP(kk) = MSE/traceP;
end
fprintf('Kalman filter \n')
fprintf('mean of MSE = %g\n',mean(MSESave))
fprintf('mean of traceP = %g\n',mean(tracePSave))
fprintf('mean of MSE/traceP = %g\n',mean(MSEdTraceP))
fprintf('std of MSE/traceP = %g\n',std(MSEdTraceP))
disp(' ')
% figure(1)
% plot([10 200], mean(MSESave)*ones(2,1),'--','Color',Color(:,2),'LineWidth',2)
% % hold on,plot([10 200], (mean(MSESave)+1*std(MSESave))*ones(2,1),'--','Color',Color(:,4),'LineWidth',2)
% % hold on, plot([10 200], (mean(MSESave)-1*std(MSESave))*ones(2,1),'--','Color',Color(:,4),'LineWidth',2)
% hold on, plot([10 200], mean(tracePSave)*ones(2,1),'--','Color',Color(:,1),'LineWidth',2)
% % hold on,plot([10 200], (mean(tracePSave)+2*std(tracePSave))*ones(2,1),'--','Color',Color(:,1),'LineWidth',2)
% % hold on, plot([10 200], (mean(tracePSave)-2*std(tracePSave))*ones(2,1),'--','Color',Color(:,1),'LineWidth',2)


%% PF
NeAll = [10 20 40 60 100 150 200];
nTries  = 50;%00; % compute MSE and traceP nTries times to average
MSEAll = zeros(length(NeAll),1);
tracePAll = zeros(length(NeAll),1);
MSEAllstd = zeros(length(NeAll),1);
tracePAllstd = zeros(length(NeAll),1);
for jj=1:length(NeAll)
    
    Ne = NeAll(jj);
    fprintf('Ne = %g\n',Ne)
    MSESave= zeros(nTries,1);
    tracePSave= zeros(nTries,1);
    MSEdTraceP = zeros(nTries,1);
    for kk=1:nTries
        % "Truth"
        xo = randn(n,1);
        xt = A*xo + Q*randn(n,1);
        y =  H*xt+sqrt(R)*randn(n,1);
        % filter
        [Xam,traceP] = myPF(y,Ne,A,Q,H,R);
        MSE = sum((Xam - xt).^2)/n;
    
        MSESave(kk) = MSE;
        tracePSave(kk) = traceP;
        MSEdTraceP(kk) = MSE/traceP;
    end
    MSEAll(jj) = mean(MSESave);
    tracePAll(jj) = mean(tracePSave);
    MSEAllstd(jj) = std(MSESave);
    tracePAllstd(jj) = std(tracePSave);
    
    MSEAll(jj)
    tracePAll(jj)
end

% save ResultsPF.mat

%%
figure(1)
hold on, plot(NeAll,MSEAll,'.-','Color',Color(:,2),'LineWidth',2,'MarkerSize',30)
hold on, plot(NeAll,tracePAll,'o-','Color',Color(:,1),'LineWidth',2,'MarkerSize',10)
box off
set(gcf,'Color','w')
set(gca,'FontSize',20)
ylabel('MSE and spread')
xlabel('N_e')






