function [xAll,yAll,MSE,traceP,z,RLocSave,RSave] = OptimizeVarPF(n,T,var_y,skip,Gap,NeAll,LocRadAll,ClocRadAll,inflAll,Save,EXIT)

%%
dt = 0.05;
t = 0:dt:T;
Steps = length(t);
%% ------------------------------------------

%% Model parameters
%% ------------------------------------------
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
muo = mean(X,2);

%% VarPF 
%% ------------------------------------------
for zz=1:length(NeAll)
    Ne = NeAll(zz);
    
    Cloc = getCov(n,ClocRadAll(zz));
    Lb = sqrtm(Cloc.*cov(X'));

    for pp = 1:length(LocRadAll)
        locRad = LocRadAll(pp);
        for tt = 1:length(inflAll)
            infl = inflAll(tt);       
            fprintf('varPF, Ne = %g, locRad = %g, infl = %g\n',Ne,locRad,infl)
            [xAll,traceP,~,RSave,RLocSave]=myVarPF(infl,locRad,Cloc,Ne,muo,Lb,z,R,H,F,Gap,Steps,dt,skip);
            MSE = mean( (xAll(:,Gap+1:Gap:end)-yAll(:,Gap+1:Gap:end)).^2 );
            if Save == 1
                Filename = strcat('varPF_Gap_',num2str(Gap),'_Ne_',num2str(Ne), ...
                    '_vary_',num2str(var_y),'_locRad_',num2str(LocRadAll(pp)),'infl',num2str(inflAll(tt)),'_skip_',num2str(skip),'.mat');
                save(Filename,'MSE','traceP','RLocSave','RSave')
            end
        end
    end
end

if EXIT ==1
    exit;
end

