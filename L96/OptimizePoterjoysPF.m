function [x,yAll,MSE,traceP,z] = OptimizePoterjoysPF(n,T,var_y,skip,Gap,NeAll,LocRadAll,aAll,Save,EXIT)
dt = 0.05;
t = 0:dt:T;
Steps = length(t);
%% ------------------------------------------

%% Model parameters
%% ------------------------------------------open 
F = 8;
%% ------------------------------------------

%% Observations
%% ------------------------------------------
H = getH(skip,n);
R = var_y*eye(size(H,1));
%% ------------------------------------------

%% Load Spin up data
%% ------------------------------------------
load(strcat('EnKFRunVary',num2str(var_y),'n',num2str(n),'Gap',num2str(Gap),'.mat'))
%% ------------------------------------------


%% Generate data
%% ------------------------------------------
yAll = model(xo,dt,Steps,F);
[z,~] = getObs(H,R,t,yAll,Gap,Steps);
%% ------------------------------------------


%% Localized OPF 
%% ------------------------------------------
for zz=1:length(NeAll)
    Ne = NeAll(zz);
    Xo = X(:,1:Ne);
    xpf = cell(Ne,1);
    for kk=1:Ne
        xpf{kk} = Xo(:,kk)';
    end
        
    for xx=1:length(LocRadAll)
        locRad = LocRadAll(xx);        
        for tt = 1:length(aAll)
            a = aAll(tt);
            fprintf('Poterjoys PF, Ne = %g, loc-rad = %g,  a = %g\n',Ne,locRad,a)
           [x,traceP] = PoterLocPF(a*Ne,locRad,Ne,xpf,z,var_y,H,F,Gap,Steps,dt,0);
           MSE = mean((x(:,Gap+1:Gap:end) - yAll(:,Gap+1:Gap:end)).^2);
           traceP = traceP'/n;
            if Save == 1
                Filename = strcat('PoterjoysPF_Gap_',num2str(Gap),'_Ne_',num2str(Ne), ...
                    '_vary_',num2str(var_y),'a',num2str(a),'LocRad',num2str(locRad),'_skip_',num2str(skip),'.mat');
                save(Filename,'MSE','traceP')
            end
        end
    end
end

if EXIT ==1
    exit;
end

