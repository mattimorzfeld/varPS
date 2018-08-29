%%
clear
close all
clc
colors


r = 0.2;
nExps = 10;
nAll = [5 10:5:50];
e = 0.1;
Ne = 1e4;

MSEwAll = zeros(length(e),1);
PwAll = zeros(length(e),1);
RAll = zeros(length(e),1);
SkewnessAll = zeros(length(e),1);

for jj=1:length(nAll)
    n = nAll(jj);
    fprintf('Dimension %g\n',n)
    MSEw_tmp = zeros(nExps,1);
    Pw_tmp = zeros(nExps,1);
    Rtmp = zeros(nExps,1);
    for kk=1:nExps  
        fprintf('     Experiment %g/%g\n',kk,nExps)
        [MSEw,Pw,R] = ExperimentUnlocalizedWeights(0*randn(n,1),r,Ne,e);
        MSEw_tmp(kk) = MSEw;
        Pw_tmp(kk) = Pw;
        Rtmp(kk) = R;
        fprintf('          R =  %g\n',R)
    end
    MSEwAll(jj) = mean(MSEw_tmp);
    PwAll(jj) = mean(Pw_tmp);
    RAll(jj) = mean(Rtmp);
end

%% 
close all
figure(1)
plot(nAll,MSEwAll,'.-','Color',Color(:,2),'LineWidth',2,'MarkerSize',30)
hold on,plot(nAll,PwAll,'o--','Color',Color(:,2),'LineWidth',2)
set(gcf,'Color','w')
set(gca,'FontSize',20)
box off


figure(2)
plot(nAll(1:length(RAll)),RAll,'.','Color',Color(:,2),'LineWidth',2,'MarkerSize',30)
set(gcf,'Color','w')
set(gca,'FontSize',20)
set(gca,'YScale','log')
xlabel('n')
ylabel('G (log-scale)')
box off

A = [nAll' ones(length(nAll),1)];
x = A\log(RAll');
figure(2)
hold on, plot(nAll,exp(x(1)*nAll+x(2)),'--','Color',Color(:,4),'LineWidth',2,'MarkerSize',25)