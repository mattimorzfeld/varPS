function out = q(x,xp,tau,y,r,e)
out = exp(-0.25/tau*(x-xp-tau*gradlogp(xp,y,r,e))^2);