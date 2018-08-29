%%
clear
close all
clc
colors

r = .2;
nExps = 5000;
n = 1;
e = 0.1:.1:1;
Ne = 40;

MSEwAll = zeros(length(e),1);
PwAll = zeros(length(e),1);
MSEnwAll = zeros(length(e),1);
PnwAll = zeros(length(e),1);
RAll = zeros(length(e),1);
SkewnessAll = zeros(length(e),1);

MSEMALAAll = zeros(length(e),1);
tracePMALAAll = zeros(length(e),1);

MSERefAll = zeros(length(e),1);
tracePRefAll = zeros(length(e),1);

MSEEnKFAll = zeros(length(e),1);
tracePEnKFAll = zeros(length(e),1);

for jj=1:length(e)
    fprintf('e = %g\n',e(jj))
    MSEw = zeros(nExps,1);
    Pw = zeros(nExps,1);
    MSEnw = zeros(nExps,1);
    Pnw = zeros(nExps,1);
    R = zeros(nExps,1);
    Skewness = zeros(nExps,1);
    
    MSEEnKF = zeros(nExps,1);
    tracePEnKF = zeros(nExps,1);
    for kk=1:nExps    
        fprintf('     Experiment %g/%g\n',kk,nExps)
        [MSEw(kk),Pw(kk),MSEnw(kk),Pnw(kk),R(kk),Skewness(kk)] ...
            = Experiment(1*randn(n,1),r,Ne,e(jj));
        [MSEEnKF(kk),tracePEnKF(kk)]=EnKF(1*randn(n,1),40,r,e(jj));
    end
    MSEEnKF
    MSEwAll(jj) = mean(MSEw);
    PwAll(jj) = mean(Pw);
    MSEnwAll(jj) = mean(MSEnw);
    PnwAll(jj) = mean(Pnw);
    RAll(jj) = mean(R);
    SkewnessAll(jj) = mean(Skewness);
    
    MSEEnKFAll(jj) = mean(MSEEnKF);
    tracePEnKFAll(jj) = mean(tracePEnKF);
end

%%
figure
plot(e,MSEwAll,'.-','Color',Color(:,6),'LineWidth',2,'MarkerSize',30)
hold on,plot(e,PwAll,'o--','Color',Color(:,6),'LineWidth',2)

hold on,plot(e,MSEnwAll,'k.-','LineWidth',2,'MarkerSize',30)
hold on,plot(e,PnwAll,'ko--','LineWidth',2)

set(gcf,'Color','w')
set(gca,'FontSize',20)
xlabel('\beta')
ylabel('MSE and spread')
box off



