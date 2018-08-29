%%
clear
close all
clc

n = 50;

%% The system
A = speye(n);
Q = speye(n);
R = speye(n);
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
    y =  H*xt+R*randn(n,1);
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


%% EnKF
Ne = 20;
nTries  = 1000; % compute MSE and traceP nTries times to average

MSESave= zeros(nTries,1);
tracePSave= zeros(nTries,1);
MSEdTraceP = zeros(nTries,1);
for kk=1:nTries
    % "Truth"
    xo = randn(n,1);
    xt = A*xo + Q*randn(n,1);
    y =  H*xt+R*randn(n,1);
    % filter
    [Xam,~,~,traceP] = myEnKF(y,Ne,A,Q,H,R,1.1);
    MSE = sum((Xam - xt).^2)/n;
    
    MSESave(kk) = MSE;
    tracePSave(kk) = traceP;
    MSEdTraceP(kk) = MSE/traceP;
end
fprintf('EnKF, ensemble size %g \n',Ne)
fprintf('mean of MSE = %g\n',mean(MSESave))
fprintf('mean of traceP = %g\n',mean(tracePSave))
fprintf('mean of MSE/traceP = %g\n',mean(MSEdTraceP))
fprintf('std of MSE/traceP = %g\n',std(MSEdTraceP))
disp(' ')


%% PF
Ne = 20;
nTries  = 1000; % compute MSE and traceP nTries times to average

MSESave= zeros(nTries,1);
tracePSave= zeros(nTries,1);
MSEdTraceP = zeros(nTries,1);
for kk=1:nTries
    % "Truth"
    xo = randn(n,1);
    xt = A*xo + Q*randn(n,1);
    y =  H*xt+R*randn(n,1);
    % filter
    [Xam,traceP] = myPF(y,Ne,A,Q,H,R);
    MSE = sum((Xam - xt).^2)/n;
    
    MSESave(kk) = MSE;
    tracePSave(kk) = traceP;  
    MSEdTraceP(kk) = MSE/traceP;
end
fprintf('PF, ensemble size %g \n',Ne)
fprintf('mean of MSE = %g\n',mean(MSESave))
fprintf('mean of traceP = %g\n',mean(tracePSave))
fprintf('mean of MSE/traceP = %g\n',mean(MSEdTraceP))
fprintf('std of MSE/traceP = %g\n',std(MSEdTraceP))
disp(' ')


%% optimal PF
Ne = 20;
nTries  = 1000; % compute MSE and traceP nTries times to average

MSESave= zeros(nTries,1);
tracePSave= zeros(nTries,1);
MSEdTraceP = zeros(nTries,1);
for kk=1:nTries
    % "Truth"
    xo = randn(n,1);
    xt = A*xo + Q*randn(n,1);
    y =  H*xt+R*randn(n,1);
    % filter
    [Xam,traceP] = myOPF(y,Ne,A,Q,H,R);
    MSE = sum((Xam - xt).^2)/n;
    
    MSESave(kk) = MSE;
    tracePSave(kk) = traceP;  
    MSEdTraceP(kk) = MSE/traceP;
end
fprintf('optimal PF, ensemble size %g \n',Ne)
fprintf('mean of MSE = %g\n',mean(MSESave))
fprintf('mean of traceP = %g\n',mean(tracePSave))
fprintf('mean of MSE/traceP = %g\n',mean(MSEdTraceP))
fprintf('std of MSE/traceP = %g\n',std(MSEdTraceP))
disp(' ')


%% Poterjoy's PF
Ne = 20;
nTries  = 100; % compute MSE and traceP nTries times to average
MSESave= zeros(nTries,1);
tracePSave= zeros(nTries,1);
MSEdTraceP = zeros(nTries,1);
locR = .1;
a = .9999;
kddm_flag = 0;
for kk=1:nTries
    % "Truth"
    xo = randn(n,1);
    xt = A*xo + Q*randn(n,1);
    y =  H*xt+R*randn(n,1);
    % filter
    [Xam,traceP] =  myPoterjoyPF(y,Ne,A,Q,H,R,locR,a,kddm_flag);
    MSE = sum((Xam - xt).^2)/n;
    
    MSESave(kk) = MSE;
    tracePSave(kk) = traceP;
    MSEdTraceP(kk) = MSE/traceP;
end
fprintf('Poterjoys PF, ensemble size %g \n',Ne)
fprintf('mean of MSE = %g\n',mean(MSESave))
fprintf('mean of traceP = %g\n',mean(tracePSave))
fprintf('mean of MSE/traceP = %g\n',mean(MSEdTraceP))
fprintf('std of MSE/traceP = %g\n',std(MSEdTraceP))

%% Poterjoy's OPF
Ne = 20;
nTries  = 100; % compute MSE and traceP nTries times to average
MSESave= zeros(nTries,1);
tracePSave= zeros(nTries,1);
MSEdTraceP = zeros(nTries,1);
locR = .1;
a = .9999;
kddm_flag = 0;
for kk=1:nTries
    % "Truth"
    xo = randn(n,1);
    xt = A*xo + Q*randn(n,1);
    y =  H*xt+R*randn(n,1);
    % filter
    [Xam,traceP] =  myPoterjoyOPF(y,Ne,A,Q,H,R,locR,a,kddm_flag);
    MSE = sum((Xam - xt).^2)/n;
    
    MSESave(kk) = MSE;
    tracePSave(kk) = traceP;
    MSEdTraceP(kk) = MSE/traceP;
end
fprintf('Poterjoys PF, ensemble size %g \n',Ne)
fprintf('mean of MSE = %g\n',mean(MSESave))
fprintf('mean of traceP = %g\n',mean(tracePSave))
fprintf('mean of MSE/traceP = %g\n',mean(MSEdTraceP))
fprintf('std of MSE/traceP = %g\n',std(MSEdTraceP))







