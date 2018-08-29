function out = F(x,y,r,e)
Mx =M(x,e);
out = .5*x.^2+(1/2/r)*(y*ones(size(Mx))-Mx).^2;