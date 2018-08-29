function [Xam,Xa,X,traceP] = myEnKF(y,Ne,A,Q,H,R,infl)
n = size(A,1);
k = size(H,1);

% initial ensemble
X = randn(n,Ne);
% localization
CL = speye(n);%gen_be_periodic(locR,1,n,1);%GetLocMatrix2(n,locR);

% forecast
X = A*X + Q*randn(n,Ne);
P = sparse(CL.*cov(X'));

Xm = mean(X,2);
Xpert = X - Xm*ones(1,Ne);
X = Xm*ones(1,Ne)+sqrt(infl)*Xpert;

K = P*H'*((H*P*H'+R)\speye(k));
Xa = X+K*(y*ones(1,Ne)+sqrt(R)*randn(k,Ne)-H*X);
Xam = Xm + K*(y-H*Xm);
traceP = trace(cov(Xa'))/n;
