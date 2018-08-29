%%
clear
close all
clc
colors

r = 1;
e = 1;

xt = 0*randn;
y = M(xt,e)+0*sqrt(r)*randn;


%% various densities
Ne=1e6;
X0 = randn(1,Ne);
X1 = M(X0,e);
w = .5*(y-X1).^2/r;
w = normalizeWeights(w);
G = mean(w.^2)/mean(w)^2;
X0rs = resampling(w,X0,Ne,1);
X1rs = M(X0rs,e);

fprintf('Posterior mean at 0: %g\n',mean(X0rs))
fprintf('Posterior cov. at 0: %g\n',cov(X0rs))
fprintf('Posterior mean at 1: %g\n',mean(X1rs))
fprintf('Posterior cov. at 1: %g\n',cov(X1rs))

figure(1)
[x, bins]=whist(X0rs,ones(Ne,1)/Ne,20);
plot(x,bins,'Color',Color(:,2),'LineWidth',2)
title('p(x_0|y)')
set(gcf,'Color','w')
set(gca,'FontSize',20)
box off

figure(2)
[x, bins]=whist(X1rs,ones(Ne,1)/Ne,50);
plot(x,bins,'Color',Color(:,2),'LineWidth',2)
title('p(x_1|y)')
set(gcf,'Color','w')
set(gca,'FontSize',20)
box off

figure(3)
[x, bins]=whist(X1,ones(Ne,1)/Ne,50);
hold on,plot(x,bins,'Color',Color(:,2),'LineWidth',2)
title('p(x_1)')
set(gcf,'Color','w')
set(gca,'FontSize',20)
box off

%% VarPS
Ne = 1e6;
xr = -4:.001:4;
Fr = F(xr,y,r,e);
[a,b]=min(Fr);
mu = xr(b);
Hess = 1+(1/r)*(Mpp(mu,e)*(M(mu,e)-y) + Mp(mu,e)^2);
C = 1/Hess;
X0varPS = mu+sqrt(C)*randn(Ne,1);
X1varPS = M(X0varPS,e);
w = F(X0varPS,y,r,e)-Fo(X0varPS,mu,C);
w = normalizeWeights(w);
G = mean(w.^2)/mean(w)^2;

figure(1)
[x, bins]=whist(X0varPS,ones(Ne,1)/Ne,30);
hold on,plot(x,bins,'Color',Color(:,4),'LineWidth',2)

figure(2)
[x, bins]=whist(X1varPS,ones(Ne,1)/Ne,1000);
hold on,plot(x,bins,'Color',Color(:,4),'LineWidth',2)
disp(' ')
fprintf('varPS mean at 0: %g\n',mean(X0varPS))
fprintf('varPS  cov. at 0: %g\n',cov(X0varPS))
fprintf('varPS  mean at 1: %g\n',mean(X1varPS))
fprintf('varPS  cov. at 1: %g\n',cov(X1varPS))
fprintf('varPS  G: %g\n',G)

%% EnKF
Ne = 1e6;
X = randn(1,Ne);
Xf = M(X,e);
P = cov(Xf);
Xm = mean(Xf);
K = P*((P+r)\1);
Xam = Xm + K*(y-Xm);
Xa = Xf+K*(y*ones(1,Ne)+1*sqrt(r)*randn(1,Ne)-Xf);
[x, bins]=whist(Xa,ones(Ne,1)/Ne,500);
figure(2)
hold on,plot(x,bins,'Color',Color(:,5),'LineWidth',2)
disp(' ')
fprintf('EnKF  mean at 1: %g\n',Xam)
fprintf('EnKF  cov. at 1: %g\n',cov(Xa))


