% Apply Kernal Density Distribution Mapping (adapted from code provided by Seth McGinnis)
%
% Input:
%        x <- first guess sample
%       xo <- original prior sample
%        w <- weights calculated using original prior   

function q = kddm(x,xo,w)

Ne = length(w);

% Save posterior mean and standard deviation
xma = sum(w.*xo);
xva = sqrt(sum(w.*(xo-xma).^2)*Ne/(Ne-1));

% Center and scale x and xo
x  = (x  - mean(x ))./sqrt(var(x ));
xo = (xo - mean(xo))./sqrt(var(xo));

% Step 1: use kernal smoothing function to approximate pdfs
[fxa, xda] = kernel_density(xo,w);
[fxf, xdf] = kernel_density(x,ones(1,Ne)/Ne);

% Step 2: use trapezoid rule to integrate pdf to get cdf
dx    = xdf(2)-xdf(1);
cdfxf = trapezoid(fxf,dx);
dx    = xda(2)-xda(1);
cdfxa = trapezoid(fxa,dx);

% Remove values with very small probability to avoid underflow error
ind1 = find(fxf<1e-4);
ind2 = find(fxf>1-1e-4);
ind = [ind1,ind2];
cdfxf(ind) = []; xdf(ind) = []; 
ind1 = find(fxa<1e-4);
ind2 = find(fxa>1-1e-4);
ind = [ind1,ind2];
cdfxa(ind) = []; xda(ind) = []; 

cdfxf = cdfxf/max(cdfxf);
cdfxa = cdfxa/max(cdfxa);

% Step 3: use spline interpolation to perform quantile matching
p = spline(xdf,cdfxf,x);
q = spline(cdfxa,xda,p);

% Step 4: re-center and scale based on posterior
q = ( q - mean(q) ).*xva./sqrt(var(q)) + xma;
