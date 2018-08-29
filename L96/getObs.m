function [z,tObs] = getObs(H,R,t,yAll,Gap,Steps)
k = size(H,1);
z = zeros(k,(Steps-1)/Gap);
nObs = size(z,2);

Lr = chol(R);
for kk=1:nObs
    jj = kk*Gap+1;
    z(:,kk) = H*yAll(:,jj)+Lr*randn(k,1);
end
tObs = t(Gap+1:Gap:end);