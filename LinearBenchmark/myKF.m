function [Xam,traceP] = myKF(y,A,Q,H,R)
n = size(A,1);
k = size(H,1);

% forecast
Xm = zeros(n,1);
P = A*speye(n)*A'+Q;

% [a,b]=meshgrid(1:n,1:n);
% mesh(a,b,P)
% view([-90 90])
% pause

K = P*H'*((H*P*H'+R)\speye(k));
Xam = Xm + K*(y-H*Xm);
P = (speye(n) - K*H)*P;

traceP = trace(P)/n;
