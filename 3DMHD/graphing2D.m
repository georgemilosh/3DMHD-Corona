load('MyColormaps')
k = 11; figure;
set(gcf, 'Position', [200 100 1100 470]); 
subplot(2,2,1); image(u(:,:,1,k),'CDataMapping','scaled'); 
axis xy, colormap(jet(256));
xlabel('x'); ylabel('y'); 
title(['Density at ' num2str(t(k))]);
colormap(mycmap);
freezeColors
colorbar; cbfreeze


subplot(2,2,2); 
image(sqrt(u(:,:,2,k).^2+u(:,:,3,k).^2),'CDataMapping','scaled'); 
axis xy; colormap(jet(256));
xlabel('x'); ylabel('y'); 
title(['x-velocity at ' num2str(t(k))]); 
hold on; h = streamslice(u(:,:,2,k),u(:,:,3,k)); 
set(h,'color','red');
hold off; 
colormap(mycmap1);
freezeColors
colorbar; cbfreeze


subplot(2,2,3); 
image(u(:,:,4,k)./u(:,:,1,k),'CDataMapping','scaled'); 
axis xy; colormap(jet(256));
xlabel('x'); ylabel('y'); 
title(['Temperature at ' num2str(t(k))]); 
colormap(hot);
freezeColors
colorbar; cbfreeze

subplot(2,2,4); 
image(u(:,:,4,k),'CDataMapping','scaled'); 
axis xy;colormap(jet(256)); 
xlabel('x'); ylabel('y'); 
title(['Pressure at ' num2str(t(k))]); 
colormap(mycmap);
freezeColors
colorbar; cbfreeze
