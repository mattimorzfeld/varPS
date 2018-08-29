function [x,yAll,MSE,traceP,z] = OptimizeEnKF(n,T,var_y,skip,Gap,NeAll,LocRadAll,inflAll,Save,EXIT)

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

%% Initial conditions/ spin up
%% ------------------------------------------
load(strcat('LongSim_n',num2str(n),'.mat'))
%% ------------------------------------------


%% Generate data
%% ------------------------------------------
yAll = model(yC(:,end),dt,Steps,F);
[z,~] = getObs(H,R,t,yAll,Gap,Steps);
%% ------------------------------------------


%% EnKF 
%% ------------------------------------------
for zz=1:length(NeAll)
    Ne = NeAll(zz);
    Xo = yC(:,randi([1 length(yC)],1,Ne));
    for pp = 1:length(LocRadAll)
        locRad = LocRadAll(pp);
        Cloc = getCov(n,locRad);
        for tt = 1:length(inflAll)
            infl = inflAll(tt);       
            fprintf('EnKF, Ne = %g, locRad = %g, infl = %g\n',Ne,locRad,infl)
            [x,traceP] = myEnKF(infl,Ne,Xo,Cloc,z,R,H,F,Gap,Steps,dt);
            MSE = mean( (x(:,Gap+1:Gap:end)-yAll(:,Gap+1:Gap:end)).^2 );
            traceP = traceP'/n;
            if Save == 1
                Filename = strcat('EnKF_Gap_',num2str(Gap),'_Ne_',num2str(Ne), ...
                    '_vary_',num2str(var_y),'_locRad_',num2str(LocRadAll(pp)),'infl',num2str(inflAll(tt)),'_skip_',num2str(skip),'.mat');
                save(Filename,'MSE','traceP')
            end
        end
    end
end

if EXIT ==1
    exit;
end

