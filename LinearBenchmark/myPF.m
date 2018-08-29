function [Xam,traceP] = myPF(y,Ne,A,Q,H,R)
n = size(A,1);
k = size(H,1);

% initial ensemble
X = randn(n,Ne);
X = A*X + 1.5* Q*randn(n,Ne);

% alpha = 1.5;

% Xm = mean(X,2);
% X = repmat(Xm,1,Ne) + alpha*( X - repmat(Xm,1,Ne) );

for kk=1:k
    x = X(kk,:);
    w = zeros(Ne,1);
    for ll = 1:Ne
        w(ll) = (y(kk) - x(ll))^2/2/R(kk,kk);
    end
    w = normalizeWeights(w);
    x = resampling(w,x,Ne,1);
    X(kk,:) = x;
end

Xam = mean(X,2);
traceP = trace(cov(X'))/n;

