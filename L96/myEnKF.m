function [xAll,traceP,X] = myEnKF(infl,Ne,X,Cloc,z,R,H,F,Gap,Steps,dt)
n = size(H,2);
nAssims = size(z,2);
k = size(H,1);
xAll = zeros(n,Steps-1);
D = zeros(k,Ne);
Lr = chol(R);

traceP = zeros(nAssims,1);
for kk=1:nAssims
    RunningMean = zeros(n,Gap+1);
    for ll=1:Ne
        trajectory = model(X(:,ll),dt,Gap+1,F);
        X(:,ll)=trajectory(:,end);
        RunningMean=RunningMean+trajectory;
        D(:,ll) = z(:,kk)+Lr*randn(k,1);
    end 
    % Intermediate time steps
    RunningMean = RunningMean/Ne;
    xAll(:,(kk-1)*Gap+1:kk*Gap+1)=RunningMean;
    
    % EnKF
    P = infl*(Cloc.*cov(X'));
    
    Xm = mean(X,2);
    Xpert = X - Xm*ones(1,Ne);
    X = Xm*ones(1,Ne)+sqrt(infl)*Xpert;
    
    K = P*H'*((H*P*H'+R)\eye(k));
    Xa = X+K*(D-H*X);
    Xam = Xm + K*(z(:,kk)-H*Xm);
    traceP(kk) = trace(cov(Xa'));
        
    X = Xa;
    xAll(:,kk*Gap+1)=Xam;%mean(Xa,2);%
end