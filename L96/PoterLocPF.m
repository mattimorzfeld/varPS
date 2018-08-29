function [xAll,traceP] = PoterLocPF(a,r,Ne,xpf,z,var_y,H,F,Gap,Steps,dt,kddm_flag)
n = size(H,2);
nAssims = size(z,2); % number of cycles
k = size(H,1); % number of obs
xAll = zeros(n,Steps);

traceP = zeros(nAssims,1);

for kk=1:nAssims
%     fprintf('Assim %g / %g\n',kk,nAssims)
    RunningMean = zeros(n,Gap+1);
    for ll=1:Ne
        trajectory = model(xpf{ll}',dt,Gap+1,F);
        RunningMean=RunningMean+trajectory;
        xpf{ll} = trajectory(:,end)';
    end
    RunningMean = RunningMean/Ne;
    xAll(:,(kk-1)*Gap+1:kk*Gap+1)=RunningMean;
    [~,xpf,eflag] = pf_update_MWR16_beta(xpf,n,Ne,H,z(:,kk)',r,a,1:k,var_y,kddm_flag);
    if eflag == 0
        fprintf('Error in localized PF.\n')
    end
    
    tmpX = zeros(n,Ne);
    for oo=1:Ne
        tmpX(:,oo) = xpf{oo}';
    end
    traceP(kk) = trace(cov(tmpX'));
end