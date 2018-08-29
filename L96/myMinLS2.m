function [x,resnorm,residual,exitflag,output,lambda,jacobian]=myMinLS2(xo,z,Gap,dt,F,H,R,mu,Lb)
func=@(x)funcF2(x,z,Gap,dt,F,H,R,mu,Lb);
options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
    'Diagnostics','off',...
    'Display','off');
n = size(H,2);
ub = inf*ones(n,1);
lb = -ub;
[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(func,xo,lb,ub,options);


