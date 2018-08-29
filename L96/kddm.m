% Apply Kernal Density Distribution Mapping (adapted from code provided by Seth McGinnis)
% This version calculates quantiles of prior explicitely and inverts posterior cdf numerically
%
% Input:
%        x <- first guess sample
%       xo <- original prior sample
%        w <- weights calculated using original prior   

function xa = kddm(x,xo,w)

warning off;

Ne = length(w);

% Save posterior mean and standard deviation
xma = sum(w.*xo);
xva = sqrt(sum(w.*(xo-xma).^2) ./ ( 1 - sum(w.^2) ) );

% Domain for defining posterior cdf
xmin = min(xo);
xmax = max(xo);
range = xmax-xmin;
incr = 2*range/(1000-1);
xd = [xmin-range/2:incr:xmax+range/2];

K = round(Ne/5);
qf = 0;
cdfxa = 0;

for i = 1:length(x)

  % Get kernel bandwidth
  sig_min = 0.001*sqrt(var(x));
  sig_max = 1.0*sqrt(var(x));
  dis = abs(x(i) - x);
  dis = sort(dis);
  sig = max(dis(K),sig_min);
  sig = min(sig,sig_max);

  % Get quantiles of prior using sum of Gaussian cdfs evaluated at prior points
  qf = qf + ( 1 + erf( ( x - x(i) )/sqrt(2)/sig ) )/2/Ne;

  % Get kernel bandwidth
  sig_min = 0.001*sqrt(var(xo));
  sig_max = 1.0*sqrt(var(xo));
  dis = abs(xo(i) - xo);
  dis = sort(dis);
  sig = max(dis(K),sig_min);
  sig = min(sig,sig_max);

  % Approximate posterior cdf
  cdfxa = cdfxa + w(i) * ( 1 + erf( (xd - xo(i))/sqrt(2)/sig ) )/2;

end

% Fix quantiles that fall out of range of domain
if min(qf) < min(cdfxa)
  r = [1:-1/(Ne-1):0];
  [a,b] = sort(x);
  qf(b) = qf(b) + r.*( min(cdfxa) - min(qf) );
end

if max(qf) > max(cdfxa)
  r = [0:1/(Ne-1):1];
  [a,b] = sort(x);
  qf(b) = qf(b) + r.*( max(cdfxa) - max(qf) );
end

% Invert cdfxa to get values at quantiles
xa = lininterp(cdfxa,xd,qf);

% Scale to fit IS mean and variance
xa = ( xa - mean(xa) ).*xva./sqrt(var(xa)) + xma;
