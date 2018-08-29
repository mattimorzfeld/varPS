function out = l(tt,rr,r,skip)


L = 2*r;
loc = skip*tt-1;
z = abs(loc-rr)/L;

if loc == 1 && rr == 40
    z = 1/L;
elseif loc == 1 && rr == 39
    z = 2/L;
elseif loc == 40 && rr == 1
    z = 1/L;
elseif loc == 40 && rr == 2
    z = 2/L;
end

if z <= 1
    out = 1 - 5/3*z^2 + 5/8*z^3 + 0.5*z^4 - 1/4*z^5;
elseif z >1 && z <= 2
    out = -2/3/z + 4 - 5*z + 5/3*z^2 + 5/8*z^3 - 1/2*z^4 + 1/12*z^5;
elseif z == 0
    out = 1;
else
    out = 0;
end