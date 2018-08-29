function H = getH(q,n)
% observe every qth variable
H = zeros(n/q,n);
k = size(H,1);
for kk=1:k
    jj = (kk-1)*q+1;
    H(kk,jj) = 1;
end