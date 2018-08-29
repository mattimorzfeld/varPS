% Function for finding beta based on effective ensemble size

function beta = find_beta_Neff(y,x,var_y,Neff)

% Original weights
beta_fg = var(y-x)/var_y;
w = exp(-(y-x).^2/2/var_y);
ws = sum(w);
w = w ./ ws;

% Values needed for calculations
Ne        = length(x);
Neff_init = 1/sum(w.^2);

if Neff == 1, beta = 1; return, end
if Neff_init < Neff || ws == 0

  % Coefficients needed for calculations
  ke = max(beta_fg,2);
  ks = 1;

  % Apply bisection method to find k
  tol  = 0.0001;

  for i = 1:100

    % Mid point
    km = (ke + ks) / 2;

    % Evaluate function at end points
    w = exp(-(y-x).^2/2/var_y/ks );
    fks = Neff - sum(w).^2 / sum(w.^2);
    if sum(w.^2) == 0, fks = 1; end
    w = exp(-(y-x).^2/2/var_y/ke );
    fke = Neff - sum(w).^2 / sum(w.^2);

    % Increase values at end point if too small
    if i == 1

      if fks < 0, break, end

      while fke*fks > 0
        ke = ke*10;
        km  = (ke + ks) / 2;
        w = exp(-(y-x).^2/2/var_y/ke);
        fke = Neff - sum(w).^2 / sum(w.^2);
        if ke > 1e5, break, end
      end
    end

    if isnan(fke), km = ke; break, end

    % Evaluate function at middle point
    w = exp(-(y-x).^2/2/var_y/km );
    fkm = Neff - sum(w).^2 / sum(w.^2);

%    disp(['iteration: ',num2str(i),', km: ',num2str(km),', fks: ',num2str(fks),', fke: ',num2str(fke)])

    % Exit criteria
    if (ke-ks)/2 < tol, break, end

    if fkm*fks > 0
      ks = km;
    else
      ke = km;
    end

  end

  % Get beta from k
  beta = km;

  w = exp(-(y-x).^2/2/var_y/beta );
  w = w ./ sum(w);

  % In extreme cases, numerical errors can lead to the wrong result
  Nf = 1./sum(w.^2);
  if abs(Neff - Nf) > 1 || isnan(Nf), beta = 99999; end
  if beta < 1, beta = 99999; end

%  Nf,Neff

%  if beta == 99999
  if beta == 0

    close all; figure(1)
    k_range = [1:0.1:1000];
    i = 0;
    for k = k_range
      i = i + 1;
      w = exp(-(y-x).^2/2/var_y/k );
      fk(i) = Neff - sum(w).^2 / sum(w.^2);
    end      
    plot(k_range,fk)

w = exp(-(y-x).^2/2/var_y );
ws = sum(w);
w
ks, ke

%    clear all; return

  end

else

  beta = 1;

end

