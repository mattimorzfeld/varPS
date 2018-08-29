function [Xam,traceP] = myOPF(y,Ne,A,Q,H,R)
n = size(A,1);
k = size(H,1);

% initial ensemble
X = randn(n,Ne);
% optimal ensemble
fX = A*X;
K = Q*H'*((H*Q*H'+R)\speye(k));
WM = R+H*Q*H';
Sigma = (speye(n)-K*H)*Q;
LSigma = chol(Sigma)';
mu = fX+K*(y*ones(1,Ne)-H*fX);
X = mu + LSigma*randn(n,Ne);

for kk=1:k
    w = zeros(Ne,1);
    for ll = 1:Ne
        w(ll) = (y(kk) - fX(kk,ll))^2/2/WM(kk,kk);
    end
    w = normalizeWeights(w);
    X(kk,:)  = resampling(w,X(kk,:) ,Ne,1);
end

Xam = mean(X,2);
traceP = trace(cov(X'))/n;

