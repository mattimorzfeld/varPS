function RunMePoterjoyPF(r,aAll,nTries)

n = 100;

%% The system
A = speye(n);
Q = speye(n);
R = r*speye(n);
H = speye(n);
k = size(H,1); % number of obs

MSESave= zeros(nTries,1);
tracePSave= zeros(nTries,1);
MSEdTraceP = zeros(nTries,1);

for ww=1:length(aAll)
    a = aAll(ww);
    
    %% Poterjoys PF
    NeAll = [10 20 40 60 100 150 200];
    MSEAll = zeros(length(NeAll),1);
    tracePAll = zeros(length(NeAll),1);
    MSEAllstd = zeros(length(NeAll),1);
    tracePAllstd = zeros(length(NeAll),1);
       
    Filename = strcat('PoterjoysPF_KDDM_r',num2str(r),'a',num2str(a),'.mat')
    for jj=1:length(NeAll)
        Ne = NeAll(jj);
        fprintf('Ne = %g\n',Ne)
        MSESave= zeros(nTries,1);
        tracePSave= zeros(nTries,1);
        MSEdTraceP = zeros(nTries,1);
        for kk=1:nTries
            % "Truth"
            xo = randn(n,1);
            xt = A*xo + Q*randn(n,1);
            y =  H*xt+sqrt(R)*randn(n,1);
            % filter
            [Xam,traceP] =  myPoterjoyPF(y,Ne,A,Q,H,R,0.01,a,1);
            MSE = sum((Xam - xt).^2)/n;
            MSESave(kk) = MSE;
            tracePSave(kk) = traceP;
            MSEdTraceP(kk) = MSE/traceP;
        end
        MSEAll(jj) = mean(MSESave);
        tmp = isnan(tracePSave);
        tmp = find(tmp==0);
        tracePSave = tracePSave(tmp);
        tracePAll(jj) = mean(tracePSave);
        MSEAllstd(jj) = std(MSESave);
        tracePAllstd(jj) = std(tracePSave);
        MSEAll(jj)
        tracePAll(jj)
    end
    save(Filename)
end







