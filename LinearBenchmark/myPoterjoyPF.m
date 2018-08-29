function [Xam,traceP] = myPoterjoyPF(y,Ne,A,Q,H,R,locR,a,kddm_flag)
n = size(A,1);
k = size(H,1);

% initial ensemble
X = randn(n,Ne);
X = A*X + Q*randn(n,Ne);
xpf = cell(Ne,1);
for ll=1:Ne
    xpf{ll} = X(:,ll)';
end

var_y = R(1,1);
[Xam,xpf,~] = pf_update_MWR16_beta(xpf,n,Ne,H,y',locR,a,1:k,var_y,kddm_flag);
tmpX = zeros(n,Ne);
for oo=1:Ne
    tmpX(:,oo) = xpf{oo}';
end
traceP = trace(cov(tmpX'))/n;
Xam = Xam';
