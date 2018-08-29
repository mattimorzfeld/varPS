% This function generates an m x m covariance matrix based on a Guassian correlation model

% Input:   
%       sigma <- standard deviation for the Gaussian correlation model
%         amp <- variance at each grid point (assumed constant at all points)
%           m <- length of the domain

function B = gen_be(sigma,amp,m,dx);

cx = m/2 + 1;
x = [1:cx,cx-1:-1:1,2:cx-1,cx:-1:1];
wlc = exp( - (x-1).*(x-1)*dx*dx/2/sigma/sigma );
for i = 1:m
  B(i,:) = wlc(m-i+2:2*m-i+1);
end


