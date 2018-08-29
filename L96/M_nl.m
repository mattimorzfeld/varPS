% This function integrates the Lorenz-96 equations forward in time 
% using the fourth order Runge-Kutta scheme
%
% Input:
%        x <- vector of size n
%       dt <- time step
%        T <- integration time
%        F <- forcing term

function x = M_nl(x,dt,T,F)

  N = length(x(1,:));

  % Create function for dx/dt
  function x1 = dxt(x0,F)

    % Create buffer zones on x for periodic domain
    x0 = [x0(1,N-1), x0(1,N), x0(1,:), x0(1,1)];

    % Place variables in vectors
    y1 = x0(1,4:end);   y2 = x0(1,3:end-1);
    y3 = x0(1,2:end-2); y4 = x0(1,1:end-3);
    x1 = (y1 - y4).*y3 - y2 + F;

  end

  for t = 1:T
    % Find RK coefficients
    k1 = dt*dxt( x(t,:)       , F );
    k2 = dt*dxt( x(t,:) + k1/2, F );
    k3 = dt*dxt( x(t,:) + k2/2, F );
    k4 = dt*dxt( x(t,:) + k3  , F );

    % Update variables
    x0 = [x(t,N-1), x(t,N), x(t,:), x(t,1)];
    y1 = x0(1,4:end);   y2 = x0(1,3:end-1);
    y3 = x0(1,2:end-2); y4 = x0(1,1:end-3);
    x(t+1,:) = y2 + k1/6 + k2/3 + k3/3 + k4/6;
  end

end
