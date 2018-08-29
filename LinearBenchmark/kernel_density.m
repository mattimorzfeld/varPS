% Function for generating a univariate kernel estimate of a pdf
%
% Input:
%       xm <- sample
%        w <- sample weights used to define posterior pdf

function [fx,x] = kernel_density(xm,w)

% Domain for perfoming the mapping
x = [min(xm)-4 : (max(xm)-min(xm)+8)/800 : max(xm)+4];

% Bandwidth for Gaussian kernels is made somewhat adaptive by
% defining it as a function of distance between variables.
% This is done by setting the bandwith equal to the distance
% between the current variable and the Kth nearest variable.

% Number of points to consider for variable bw
% (5 seems to be a good choice, but other values can be used)
K = floor(length(xm)/5);

fx = 0;
for i = 1:length(xm)

  % Distance between current data point and all other points
  dis = abs(xm(i) - xm);  

  % Sort distances
  dis = sort(dis);

  % Use distance between Kth nearest point and current point for bw
  sig = max(dis(K),0.01);

  % Sum Gaussian kernels
  fx = fx + w(i)*exp(- 0.5 * (x-xm(i)).^2/sig^2 )/sqrt(2*pi)/sig;

end
