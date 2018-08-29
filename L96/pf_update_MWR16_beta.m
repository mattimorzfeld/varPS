% Function for performing the particle filter update step

function [xmpf,xa,eflag] = pf_update_MWR16_beta(xpf,Nx,Ne,H,Y,roi,Neff,oborder,var_y,kddm_flag)

% Generate correlation matrix for localization
C = gen_be_periodic(roi,1,Nx,1);

% Save prior ensemble
xpfo = xpf;

% Initialize weighting matrix with ones
wo = ones(Ne,Nx)/Ne;

for i = oborder % Observation loop

  % Set operator for observation i
  Hi = H(i,:);

  % Observation space priors 
  clear hx hxo
  for n = 1:Ne
    hxo(n) = Hi*xpfo{n}';
    hx(n)  = Hi*xpf{n}';
  end

  beta = find_beta_Neff(Y(i),hxo,var_y,Neff);
  
  % Set vector of localization coefficients
  loc = Hi*C;

  % Update weights: w is the set of scalar weights calculated using
  % the current observation and particles that have been updated
  % by all previous observations; wo is the set of vector weights
  % calculated recursively using all observations up to the current
  % one, and a prior ensemble that has not been updated by previous
  % observations in the data assimilation cycle.

  % Calculate w
  d1 = (Y(i)-hx)/sqrt(2*var_y*beta);   % Normalized innovation
  w = exp( - d1.*d1 )./sqrt(2*pi);        % p(y|x_n)

  % Calculate wo
  d2 = (Y(i)-hxo)/sqrt(2*var_y*beta);   % Normalized innovation
  wn = exp( - d2.*d2 )./sqrt(2*pi);         % p(y|x_n)

  % Keep normalization
  ws = sum(w.*(Hi*wo'));

  % Ignore ob if weight collapse still occurs 
  wt = [wo(:,:)*Hi']'.*( ( wn - 1 ) + 1  );
  if sum(wt) == 0
      xmpf = 0;
      for n = 1:Ne
          xmpf = xmpf+xpfo{n};
      end
      xmpf = xmpf/Ne;
      for n = 1:Ne
          xa{n} = xpfo{n};
      end
      eflag = 1;
  else
      % Set vector of localization coefficients and localize weights
      loc = Hi*C;
      for n = 1:Ne
          wo(n,:) = wo(n,:).*( ( wn(n) - 1 ).*loc + 1  );
      end
      
      % Normalize weights
      w   = w/sum(w);
      wos = sum(wo);
      
      for n = 1:Ne
          wo(n,:) = wo(n,:)./wos;
      end
      % NOTE: For high-dimensional systems, it helps to normalize wo after each
      %       observation is assimilated. Underflow error can sometimes occur with
      %       densely-spaced observations.
      
      % Calculate posterior mean using vector weights and original prior particles
      xmpf = 0;
      for n = 1:Ne
          xmpf = xmpf + wo(n,:).*xpfo{n};
      end
      
      % Calculate posterior variance using vector weights and original prior particles
      var_a = 0;
      for n = 1:Ne
          var_a = var_a + wo(n,:).*(xpfo{n}-xmpf).*(xpfo{n}-xmpf);
      end
      
      for j = 1:Nx
          var_a(j) = var_a(j) ./ ( 1 - sum(wo(:,j).^2) );
      end
      if  sum(isnan(var_a)) ~= 0
          xmpf = 0;
          for n = 1:Ne
              xmpf = xmpf+xpfo{n};
          end
          xmpf = xmpf/Ne;
          for n = 1:Ne
              xa{n} = xpfo{n};
          end
          eflag = 1;      
      else
          % Error check
          if sum(sum(isnan(xmpf))) > 0
              eflag = 0; return
          else
              eflag = 1;
          end
          
          % Apply deterministic resampling
          ind = sampling(hx,w',Ne);
          
          % Calculate c term for weight update equation
          alpha = 0.99;
          c = Ne*(1-loc*alpha)./loc./ws./alpha;
          
          % Calculate r1 and r2 coefficients for weight update equation
          r1 = 0; r2 = 0;
          
          for n = 1:Ne
              r1 = r1 + (  xpf{ind(n)}-xmpf + c .* ( xpf{n}-xmpf ) ).^2;
              r2 = r2 + ( (xpf{ind(n)}-xmpf)./c +  ( xpf{n}-xmpf ) ).^2;
          end
          r1 = sqrt((Ne-1)*var_a./r1);
          if sum(isnan(r1))~=0
              r1 = zeros(1,Nx);
          end
          r2 = sqrt((Ne-1)*var_a./r2);
          if sum(isnan(r2))~=0
              r2 = zeros(1,Nx);
          end
          
          % Generate localized posterior particles
          for n = 1:Ne
              xa{n} = xmpf + r1.*(xpf{ind(n)} - xmpf) + r2.*(xpf{n} - xmpf);
          end
          
          % Adjust posterior mean and variance to correct for sampling errors
          vs = 0; pfm = 0;
          for n = 1:Ne
              pfm = pfm + xa{n}./Ne;
          end
          
          for n = 1:Ne
              vs = vs + (xa{n} - pfm).*(xa{n} - pfm)./(Ne-1);
          end
          
          correction = sqrt(var_a./vs);
          %       fprintf('correction = %g\n',correction)
          if sum(isnan(correction))~=0
              correction = zeros(1,Nx);
          end
          for n = 1:Ne
              xa{n} = xmpf + (xa{n} - pfm).*correction;
          end
          
          % Make posterior ensemble the new prior ensemble
          xpf = xa;
      end
  end
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
    q = kddm(x,xo,wo(:,j)');
    for n = 1:Ne
      xa{n}(j) = q(n);
    end
  end
end
