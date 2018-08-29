function myerrorCloud(x,y,z,Color1,Color2)
h = patch([z;z(end);z(end:-1:1);z(1)],[x+y;x(end)-y(end);x(end:-1:1)-y(end:-1:1);x(1)+y(2)],'g');
set(h,'FaceColor',Color1)
set(h,'EdgeColor',Color1)
% set(h,'FaceAlpha',0.4)
% set(h,'EdgeAlpha',0.4)
hold on, plot(z,x,'Color',Color2,'LineWidth',2)