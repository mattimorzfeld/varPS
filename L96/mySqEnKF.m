function [xAll,traceP,X] = mySqEnKF(infl,Ne,X,Cloc,z,R,H,F,Gap,Steps,dt)
n = size(H,2);
nAssims = size(z,2);
k = size(H,1);
xAll = zeros(n,Steps-1);
D = zeros(k,Ne);
Lr = chol(R);

traceP = zeros(nAssims,1);
for kk=1:nAssims
%     fprintf('Assim %g / %g\n',kk,nAssims)
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

    P = infl*(Cloc.*cov(X'));
    
    Xm = mean(X,2);
    Xpert = X - Xm*ones(1,Ne);
    X = Xm*ones(1,Ne)+sqrt(infl)*Xpert;
    
    K = P*H'*((H*P*H'+R)\eye(k));
    Xam = Xm + K*(z(:,kk)-H*Xm);
    xAll(:,kk*Gap+1)=Xam;
   
    Z = (sqrt(infl)/sqrt(Ne-1))*Xpert;

    tmp = Z'*H'*(R\(H*Z));
    tmp = .5*(tmp+tmp');
    [E,OM] = eig(tmp);
    T = E*(sqrtm(eye(Ne)+OM)\E');
    Za = Z*T;
    X = repmat(Xam,1,Ne)+sqrt(Ne-1)*Za;
    
    traceP(kk) = trace(cov(X'));
end