function CL = GetLocMatrix2( Nx,dsm )

%% spatial domain
xS = -pi;
xF = pi;
dx = (xF - xS)/(Nx - 0);
x = [xS:dx:xF-dx]';

%% create sinusoid matrix (basis functions)
% ... because we use sinusoids our correlation matrix will be periodic!

E = zeros(Nx);
zj = (2*pi/Nx:2*pi/Nx:2*pi)';
e0 = ones(Nx,1)./sqrt(Nx);
E(:,1) = e0;
for i = 1:Nx/2-1,
  tt = cos(i*zj);
  tn = sqrt(tt'*tt);
  E(:,2*i) = tt./tn;
  tt = sin(i*zj);
  tn = sqrt(tt'*tt);
  E(:,2*i+1) = tt./tn;
end
tt = cos(Nx*zj/2);
tn = sqrt(tt'*tt);
E(:,Nx) = tt./tn;

%% create eigenvalues

g = ones(1,Nx);
b = 1./dsm^2;
for i = 1:(Nx-2)/2
  g(2*i)   = exp(-b*i^2);
  g(2*i+1) = g(2*i);
end
g(Nx) = exp(-b*(Nx/2)^2);

a = Nx/sum(g);
g = a*g;
G = diag(g);

CL = E*G*E';
CL = CL/max(max(CL));

end