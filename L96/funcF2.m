function f = funcF2(x,z,Gap,dt,F,H,R,mu,Lb)
n = length(x);
X= mu+Lb*x;
trajectory = model(X,dt,Gap+1,F);
X = trajectory(:,end);
k = size(H,1); % number of obs
f = zeros(k,1);
for kk=1:k
    f(kk) = (H(kk,:)*X-z(kk))/sqrt(2*R(kk,kk));
end
f = real([f;sqrt(.5)*x]);
% disp('Lb = '), Lb
% disp('mu = '),mu'
% disp('f = '),f'