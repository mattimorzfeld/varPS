% Function that performs pf sampling step
%
% Input:
%        w <- weights used for determining sampling
%       Ne <- ensemble size

function ind = sampling(x,w,Ne);

% Sort sample
[a,b] = sort(x);

% Apply deterministic sampling by taking value at every 1/Ne quantile
cum_weight = [0;cumsum(w(b))];
base = 1/Ne/2;
clear ind
ind  = [];
for n = 1:Ne
  frac = base + (n - 1)/Ne;
  for k = 2:Ne+1
    if ( cum_weight(k-1) < frac ) && ( frac <= cum_weight(k) )
      ind(n) = k-1;
      continue
    end
  end
end

% Unsort indices
ind = b(ind);

% Replace removed particles with duplicated particles
ind2 = ind*0;
for n = 1:Ne
  if sum(ind==n) ~= 0
    ind2(n) = n;
    dum = find(ind==n);
    ind(dum(1)) = [];
  end
end

ind0 = find(ind2==0);
ind2(ind0) = ind;
ind = ind2;
if isempty(ind)
    ind = 1:Ne;
end
