function out = normpdf(x,xm,s)
out =exp(- (x-xm).^2/2/s^2 )/s;
% out = sqrt(2*pi*s)\exp(- (y-x)^2/2/s );