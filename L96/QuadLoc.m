function f = QuadLoc(ind,n,locR)
d1 = abs(ind-(1:n));
d2 = abs(ind+n -(1:n));
d3 = abs(ind-n -(1:n));

f1 = exp(-(d1/locR).^2);
f2 =  exp(-(d2/locR).^2);
f3 =  exp(-(d3/locR).^2);

f = max([f1;f2;f3]);