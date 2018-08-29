function myerrorbar(x,y,z,Color,DX,Marker,MarkerSize)
X = x;
Y = y;
Z = z;
for kk=1:length(X)
    x=X(kk); y = Y(kk); z = Z(kk);
    hold on, plot(z,x,Marker,'Color',Color,'MarkerSize',MarkerSize,'LineWidth',2)
    hold on, plot([z z],[x+y,x-y],'-','Color',Color,'LineWidth',2)
    hold on, plot([z+DX z-DX],[x+y,x+y],'-','Color',Color,'LineWidth',2)
    hold on, plot([z+DX  z-DX],[x-y,x-y],'-','Color',Color,'LineWidth',2)
end