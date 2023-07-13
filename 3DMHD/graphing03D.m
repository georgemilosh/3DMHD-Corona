% time to be drawn

figure
set(gcf, 'Position', [0 100 1366 450]);
subplot(1,2,1)
hsurfaces = slice(x,y,z,u0(:,:,:,1),round(size(x,2)/2),...
    round(size(y,1)/2),1);%round(size(z,3)/2)
set(hsurfaces,'EdgeColor','none')
colormap(jet(256)); 
title('N at time = 0')
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar;
view(0,0)
axis equal


subplot(1,2,2)
hsurfaces = slice(x,y,z,u0(:,:,:,5),round(size(x,2)/2),...
    round(size(y,1)/2),1);%round(size(z,3)/2)
set(hsurfaces,'EdgeColor','none')
colormap(jet(256)); 
title('P at time = 0')
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar;
view(0,0)
axis equal


figure 
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
set(gcf, 'Position', [0 100 1366 350]);

subplot(1,3,1)

hsurfaces = slice(x,y,z,u0(:,:,:,2),round(size(x,2)/2),...
    round(size(y,1)/2),1);%round(size(z,3)/2)
set(hsurfaces,'EdgeColor','none')
colormap(jet(256)); 
title('V_x at time = 0')
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar;
view(0,0)
axis equal


subplot(1,3,2)
hsurfaces = slice(x,y,z,u0(:,:,:,3),round(size(x,2)/2),...
    round(size(y,1)/2),1);%round(size(z,3)/2)
set(hsurfaces,'EdgeColor','none')
colormap(jet(256)); 
title('V_y at time = 0')
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar;
view(0,0)
axis equal

subplot(1,3,3)
hsurfaces = slice(x,y,z,u0(:,:,:,4),round(size(x,2)/2),...
    round(size(y,1)/2),1);%round(size(z,3)/2)
set(hsurfaces,'EdgeColor','none')
colormap(jet(256)); 
title('V_z at time = 0')
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar;
view(0,0)
axis equal


figure 
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
set(gcf, 'Position', [0 100 1366 350]);

subplot(1,3,1)

hsurfaces = slice(x,y,z,u0(:,:,:,6),round(size(x,2)/2),...
    round(size(y,1)/2),1);%round(size(z,3)/2)
set(hsurfaces,'EdgeColor','none')
colormap(jet(256)); 
title('B_x at time = 0')
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar;
view(0,0)
axis equal


subplot(1,3,2)
hsurfaces = slice(x,y,z,u0(:,:,:,7),round(size(x,2)/2),...
    round(size(y,1)/2),1);%round(size(z,3)/2)
set(hsurfaces,'EdgeColor','none')
colormap(jet(256)); 
title('B_y at time = 0')
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar;
view(0,0)
axis equal

subplot(1,3,3)
hsurfaces = slice(x,y,z,u0(:,:,:,8),round(size(x,2)/2),...
    round(size(y,1)/2),1);%round(size(z,3)/2)
set(hsurfaces,'EdgeColor','none')
colormap(jet(256)); 
title('B_z at time = 0')
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar;
view(0,0)
axis equal
%saveas(gcf, 'myfilename.jpg');
