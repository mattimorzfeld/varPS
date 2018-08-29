% Function for performing the particle filter update step

function [xmpf,xa,eflag] = opf_update(xpf,xpfo,fXnm1,Nx,Ne,H,Y,roi,alpha,oborder,var_y,kddm_flag,Q)

% Generate correlation matrix for localization
C = gen_be_periodic(roi,1,Nx,1);

% Save prior ensemble
% xpfo = xpf;

% Initialize weighting matrix with ones
wo = ones(Ne,Nx);

WM = var_y*eye(length(oborder))+H*Q*H';

for i = oborder % Observation loop

  % Set operator for observation i
  Hi = H(i,:);
  
  % Observation space priors 
  clear hx hxo
  for n = 1:Ne
    hxo(n) = Hi*xpfo{n}';
    hx(n)  = Hi*xpf{n}';
  end

  % Set vector of localization coefficients and multiply by alpha
  loc = Hi*C.*alpha;

  % Update weights: w is the set of scalar weights calculated using
  % the current observation and particles that have been updated
  % by all previous observations; wo is the set of vector weights
  % calculated recursively using all observations up to the current
  % one, and a prior ensemble that has not been updated by previous
  % observations in the data assimilation cycle.

  clear w
  for n = 1:Ne

    % Calculate w
    wn = exp(-.5*(Y(i) - Hi*fXnm1(:,n))^2/ WM(i,i));
    w(n) = ( wn - 1 )*alpha + 1;

    % Calculate wo
    wn = exp(-.5*(Y(i) - Hi*fXnm1(:,n))^2/ WM(i,i));
    wo(n,:) = wo(n,:).*( ( wn - 1 ).*loc + 1 );

    % NOTE: For Gaussian observation errors, exp( (y-hx)/sqrt(2*var_y) )
    %       should be divided by sqrt(2*pi*var_y). You may notice that
    %       var_y is left out of the denominator for wn in the weight
    %       calculations. This factor does not matter for the standard
    %       bootstrap filter, since the weights are normalized by their
    %       sum. But for the current filter, using a smaller or larger
    %       observation error variance will change the effects of the
    %       localization. This dependence is removed by ignoring var_y
    %       in the denominator.

  end

  % Normalize weights
  ws  = sum(w);
  w   = w/ws;
  wos = sum(wo);

  % NOTE: For high-dimensional systems, it helps to normalize wo after each
  %       observation is assimilated. Underflow error can sometimes occur with
  %       densely-spaced observations.

  % Calculate posterior mean using vector weights and original prior particles
  xmpf = 0;
  for n = 1:Ne
    xmpf  = xmpf + (wo(n,:)./wos).*xpfo{n};
  end

  % Calculate posterior variance using vector weights and original prior particles
  var_a = 0;
  for n = 1:Ne
    var_a = var_a + (wo(n,:)./wos).*(xpfo{n}-xmpf).*(xpfo{n}-xmpf)*Ne/(Ne-1);
  end
    
   if sum(sum(isnan(xmpf))) > 0
       fprintf('isnan detected\n')
       sum(isnan(xmpf))
   end
  % Error check
  if sum(sum(isnan(xmpf))) > 0
    eflag = 0; return
  else
    eflag = 1;
  end

  % Apply deterministic resampling
  ind = sampling(hx,w',Ne);

  % Calculate c term for weight update equation
  c = Ne*(1-loc)./loc./ws;

  % Calculate r1 and r2 coefficients for weight update equation
  r1 = 0; r2 = 0;
  for n = 1:Ne
    r1 = r1 + ( xpf{ind(n)}-xmpf + c .* ( xpf{n}-xmpf ) ).^2;
    r2 = r2 + ( (xpf{ind(n)}-xmpf)./c + ( xpf{n}-xmpf ) ).^2;
  end
  r1 = sqrt((Ne-1)*var_a./r1);
  r2 = sqrt((Ne-1)*var_a./r2);

  % Generate localized posterior particles
  for n = 1:Ne
    xa{n} = xmpf + r1.*(xpf{ind(n)} - xmpf) + r2.*(xpf{n} - xmpf);
  end

  % Adjust posterior mean and variance to correct for sampling errors
  vs = 0; pfm = 0; var_p = 0; vm = 0; pm = 0;
  for n = 1:Ne
    pfm = pfm + xa{n}./Ne;
     vm = vm  + xpfo{n}./Ne;
     pm = pm  + xpf{n}./Ne;
  end
  for n = 1:Ne
    var_p = var_p + (xpfo{n} - vm).*(xpfo{n} - vm)./(Ne-1);
    % From updated members
    vs = vs + (xa{n} - pfm).*(xa{n} - pfm)./(Ne-1);
  end

  correction = sqrt(var_a)./sqrt(vs);
  for n = 1:Ne
    xa{n} = xmpf + (xa{n} - pfm).*correction;
  end

  % Make posterior ensemble the new prior ensemble
  xpf = xa;

end

% Use KDDM to update posterior samples
if kddm_flag
  for j = 1:Nx % State variable loop
    clear x xo
    % Place prior and current posterior representation of variable j in vectors.
    for n = 1:Ne
      x(n)  = xpf{n}(j);
      xo(n) = xpfo{n}(j);
    end
    q = kddm(x,xo,wo(:,j)'/wos(j));
    for n = 1:Ne
      xa{n}(j) = q(n);
    end
  end
end
