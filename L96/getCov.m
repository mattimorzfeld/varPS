function C = getCov(n,L)
% squared exponential
C = zeros(n);
for ii=1:n
    for jj=ii:n
        dist = min(abs(ii-jj), mod(-abs(ii-jj),n));
        C(ii,jj) = exp(-(dist)^2/(2*L^2));
    end
end
C = (C+C')-diag(diag(C));


%% Old version: sin based covariance
%C = getCov(n,T,L)
% for jj = 1:n
%     for ii=jj:n
%         C(jj,ii)=exp(-2*sin(abs(ii-jj)/(n-1)*2*pi/2)^2/L^2);
%     end
% end
% C = C+C'-diag(diag(C));