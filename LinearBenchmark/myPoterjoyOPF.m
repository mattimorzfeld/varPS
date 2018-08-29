function [Xam,traceP] = myPoterjoyOPF(y,Ne,A,Q,H,R,locR,a,kddm_flag)
n = size(A,1);
k = size(H,1);

% initial ensemble
X = randn(n,Ne);
% get optimal ensemble
fX = A*X;
% model forecast
Xf = fX +Q*randn(n,Ne);
% Gain
K = Q*H'*((H*Q*H'+R)\eye(k));
% WM = R+H*Q*H';
Sigma = (eye(n)-K*H)*Q;
LSigma = chol(Sigma)';

mu = fX+K*(y*ones(1,Ne)-H*fX);
X = mu + LSigma*randn(n,Ne);

xpf = cell(Ne,1);
xpfo = cell(Ne,1);
for ll=1:Ne
    xpf{ll} = X(:,ll)';
    xpfo{ll} = Xf(:,ll)';
end

var_y = R(1,1);
[~,xpf,~] = opf_update(xpfo,xpf,fX,n,Ne,H,y',locR,a,1:k,var_y,kddm_flag,Q);
% [~,xpf,~] = pf_update(xpfo,n,Ne,H,y',locR,a,1:k,var_y,kddm_flag);


tmpX = zeros(n,Ne);
for oo=1:Ne
    tmpX(:,oo) = xpf{oo}';
end
traceP = trace(cov(tmpX'))/n;
Xam = mean(tmpX,2);
