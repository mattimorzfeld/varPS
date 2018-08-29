% Linear interpolation

function xi = lininterp(fx,x,fxi);

for n = 1:length(fxi)

  [dum,m] = min( fx < fxi(n) );

  if fxi(n) >= max(fx)
    xi(n) = max(x);
    continue
  end

  if fxi(n) <= min(fx)
    xi(n) = min(x);
    continue
  end

  if fxi(n) == fx(m) || m == 1
    xi(n) = x(m);
  else
    if fx(m) > fxi(n), m = m - 1; end
    w1 = ( fx(m+1) - fxi(n) ) / ( fx(m+1) - fx(m) );
    w2 = ( fxi(n) -  fx(m)  ) / ( fx(m+1) - fx(m) );
    xi(n) = w1*x(m) + w2*x(m+1);
  end

end
