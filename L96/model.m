function yAll = model(y,dt,Steps,F)
n = length(y);
yAll = zeros(n,Steps);
yAll(:,1) = y;
for kk=2:Steps
    k1 = dt * fL40(y, F);
    k2 = dt * fL40(y + 0.5 * k1, F);
    k3 = dt * fL40(y + 0.5 * k2, F);
    k4 = dt * fL40(y + k3, F);
    y = y + (k1 + 2 * (k2 + k3) + k4) / 6;
    yAll(:,kk) = y;
end
end