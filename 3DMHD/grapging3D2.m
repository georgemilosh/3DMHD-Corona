
% time to be drawn
time = 1;

figure
set(gcf, 'Position', [0 100 1366 500]);
subplot(1,2,1)
hsurfaces = slice(x,y,z,u(:,:,:,1),round(size(x,2)/2),...
    round(size(y,1)/2),1);%round(size(z,3)/2)
set(hsurfaces,'EdgeColor','none')
colormap(jet(256)); 
title(['N at time = ' num2str(t(time))])
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar;
view(0,0)
axis tight


subplot(1,2,2)
hsurfaces = slice(x,y,z,u(:,:,:,5),round(size(x,2)/2),...
    round(size(y,1)/2),1);%round(size(z,3)/2)
set(hsurfaces,'EdgeColor','none')
colormap(jet(256)); 
title(['P at time = ' num2str(t(time))])
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar;
view(0,0)
axis tight


figure 
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
set(gcf, 'Position', [0 100 1366 500]);

subplot(1,2,1)

hsurfaces = slice(x,y,z,sqrt(u(:,:,:,2).^2+u(:,:,:,3).^2+u(:,:,:,4).^2),...
    round(size(x,2)/2),...
    round(size(y,1)/2),1);%round(size(z,3)/2)
set(hsurfaces,'EdgeColor','none')
colormap(jet(256)); 
title(['V at time = ' num2str(t(time))])
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar;
view(0,0)
axis tight
hold on
[sx,sy,sz] = meshgrid(1:round(size(x,2)/10):round(size(x,2)/2),...
    round(size(y,1)/4):round(size(y,1)/4):(size(y,1)-round(size(y,1)/4)),1);
hlines = streamline(x,y,z,u(:,:,:,2),u(:,:,:,3),u(:,:,:,4),sx,sy,sz);
set(hlines,'LineWidth',2,'Color','r')
hold off

subplot(1,2,2)

hsurfaces = slice(x,y,z,sqrt(u(:,:,:,6).^2+u(:,:,:,7).^2+u(:,:,:,8).^2),...
    round(size(x,2)/2),...
    round(size(y,1)/2),1);%round(size(z,3)/2)
set(hsurfaces,'EdgeColor','none')
colormap(jet(256)); 
title(['B at time = ' num2str(t(time))])
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar;
view(0,0)
axis tight
hold on
[sx,sy,sz] = meshgrid(1:round(size(x,2)/10):round(size(x,2)/2),...
    round(size(y,1)/4):round(size(y,1)/4):(size(y,1)-round(size(y,1)/4)),1);
hlines = streamline(x,y,z,u(:,:,:,6),u(:,:,:,7),u(:,:,:,8),sx,sy,sz);
set(hlines,'LineWidth',2,'Color','r')
hold off


