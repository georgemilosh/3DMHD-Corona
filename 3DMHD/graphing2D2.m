load('MyColormaps')
k = 10; 

R = 10e9/size(u(:,:,1,1),2);
[x y] = meshgrid(linspace(0,R*size(u(:,:,1,1),2),size(u(:,:,1,1),2)),linspace(0,R*size(u(:,:,1,1),1),size(u(:,:,1,1),1)));
RVa = R/Va;
figure;
%set(gcf, 'Position', [50 50 850 630]); 
image(x(1,:),y(:,1),n0*u(:,:,1,k),'CDataMapping','scaled'); 
axis xy;
xlabel('x [cm]'); ylabel('y [cm]');
title(['Density [cm^-3]; t = ' num2str(round(RVa*t(k))) ' s']);
colormap(jet);
colorbar;

figure;
image(x(1,:),y(:,1),Va*sqrt(u(:,:,2,k).^2+u(:,:,3,k).^2),'CDataMapping','scaled'); 
axis xy;
caxis([0 3e7]);
xlabel('x [cm]'); ylabel('y [cm]');
title(['Velocity [cm/s]; t = ' num2str(round(RVa*t(k))) ' s']);
hold on; h = streamslice(x,y,u(:,:,2,k),u(:,:,3,k)); 
set(h,'color','red');
hold off; 
colormap(mycmap1);
colorbar;

figure;
image(x(1,:),y(:,1),T0*beta^-1*u(:,:,4,k)./u(:,:,1,k)/1.6e-12,'CDataMapping','scaled'); 
axis xy; colormap(jet(256));
xlabel('x [cm]'); ylabel('y [cm]');
title(['Temperature [eV]; t = ' num2str(round(RVa*t(k))) ' s']);
colormap(hot);
colorbar;

figure;
image(x(1,:),y(:,1),B0*sqrt(u(:,:,5,k).^2+u(:,:,6,k).^2),'CDataMapping','scaled'); 
axis xy;colormap(jet(256)); 
xlabel('x [cm]'); ylabel('y [cm]');
title(['Magnetic Field [Gauss]; t = ' num2str(round(RVa*t(k))) ' s']);
hold on; h1 = streamslice(x,y,u(:,:,5,k),u(:,:,6,k)); 
set(h1,'color','red');
hold off; 
colormap(pink);
colorbar;
