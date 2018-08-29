% Function for finding alpha based on effective ensemble size

function alpha = find_alpha_Neff(y,x,var_y,Neff)

% Original weights
w = exp(-(y-x).^2/2/var_y )/sqrt(2*pi);
ws = sum(w);
w = w ./ ws;

% Values needed for calculations
Ne        = length(x);
Neff_init = 1/sum(w.^2);

if Neff_init < Neff || ws == 0

  % Coefficients needed for calculations
  a = exp(-(y-x).^2/2/var_y)/sqrt(2*pi);
  ae = 1;
  as = 0;

  % Apply bisection method to find k
  tol  = 0.0001;

  for i = 1:100

    % Mid point
    am = (ae + as) / 2;

    % Evaluate function at three points
    fas = Neff*sum( ( (a - 1)*as + 1 ).^2 ) - ( Ne + as * sum( a - 1 ) ).^2;
    fae = Neff*sum( ( (a - 1)*ae + 1 ).^2 ) - ( Ne + ae * sum( a - 1 ) ).^2;
    fam = Neff*sum( ( (a - 1)*am + 1 ).^2 ) - ( Ne + am * sum( a - 1 ) ).^2;

    % Exit criteria
    if (ae-as)/2 < tol, break, end

    if fam*fas > 0
      as = am;
    else
      ae = am;
    end

%    disp(['iteration: ',num2str(i),' am: ',num2str(am)])

  end

  alpha = am;

  w = ( a - 1 ) * alpha + 1;
  w = w ./ sum(w);

  % Numerical errors can still lead to the wrong result
  Nf =  1./sum(w.^2);
  if abs(Neff - Nf) > 0.1 || isnan(Nf), alpha = 0; end

%  disp(['Starting Neff ',num2str(Neff_init),' Target Neff: ',num2str(Neff),' New Neff: ',num2str(1/sum(w.^2))])
%  disp(['Alpha: ',num2str(alpha)])

else

  alpha = 1;

end
