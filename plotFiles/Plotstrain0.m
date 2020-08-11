function Plotstrain0(FSubpb3,x,y,sizeOfImg,CurrentImg,OrigDICImgTransparency)
  
warning off;

M = size(x,1); N = size(x,2);
u_x = FSubpb3(1:4:end); v_x = FSubpb3(2:4:end);
u_y = FSubpb3(3:4:end); v_y = FSubpb3(4:4:end);
 
u_x = reshape(u_x,M,N); v_x = reshape(v_x,M,N);
u_y = reshape(u_y,M,N); v_y = reshape(v_y,M,N);

% imagesc([x(1,1) x(end,1)], [y(1,1) y(1,end)], flipud(g)); hold on;
if M < 9,x2 = x(:,1)'; else x2 = interp(x(:,1)',4); end
if N < 9,y2 = y(1,:); else y2 = interp(y(1,:),4); end

z_exx = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u_x,M*N,1),x2,y2);

 


%figure; 
% contourf(x2,y2,z,10,'EdgeColor','none','LineStyle','none')
% set(gca,'fontSize',18); view(2);
% set(gca,'ydir','normal');
% title('Strain $e_{xx}$','FontWeight','Normal','Interpreter','latex');
% axis tight; axis equal; colorbar; 
% if x(M,N) < 200
%     set(gca,'XTick',[]);
% end
% if y(M,N) < 200
%     set(gca,'YTick',[]);
% end
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% set(gcf,'color','w');
% a = gca; a.TickLabelInterpreter = 'latex';
% b = colorbar; b.TickLabelInterpreter = 'latex';
% 
% colormap jet
% box on

% alpha(c,.5)

figure;  
surf(x2,sizeOfImg(2)+1-y2,z_exx,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('Strain $e_{xx}$','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
% alpha(c,.5)
 



z_exy = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(0.5*(v_x+u_y),M*N,1),x2,y2);
% figure; 
% contourf(x2,y2,z,10,'EdgeColor','none','LineStyle','none')
% set(gca,'fontSize',18); view(2); 
% title('Strain $e_{xy}$','FontWeight','Normal','Interpreter','latex');
% axis tight; axis equal; colorbar; 
% if x(M,N) < 200
%     set(gca,'XTick',[]);
% end
% if y(M,N) < 200
%     set(gca,'YTick',[]);
% end
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% set(gcf,'color','w');
% a = gca; a.TickLabelInterpreter = 'latex';
% b = colorbar; b.TickLabelInterpreter = 'latex';
% colormap jet
% box on


figure;  
surf(x2,sizeOfImg(2)+1-y2,z_exy,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('Strain $e_{xy}$','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
% alpha(c,.5)


z_eyy = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v_y,M*N,1),x2,y2);
% figure;  
% contourf(x2,y2,z,10,'EdgeColor','none','LineStyle','none')
% set(gca,'fontSize',18); view(2); 
% title('Strain $e_{yy}$','FontWeight','Normal','Interpreter','latex');
% axis tight; axis equal; colorbar; 
% if x(M,N) < 200
%     set(gca,'XTick',[]);
% end
% if y(M,N) < 200
%     set(gca,'YTick',[]);
% end
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% set(gcf,'color','w');
% a = gca; a.TickLabelInterpreter = 'latex';
% b = colorbar; b.TickLabelInterpreter = 'latex';
% colormap jet
% box on


figure;  
surf(x2,sizeOfImg(2)+1-y2,z_eyy,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('Strain $e_{yy}$','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
% alpha(c,.5)



figure;  
z_maxshear = sqrt((0.5*(z_exx-z_eyy)).^2 + z_exy.^2);
surf(x2,sizeOfImg(2)+1-y2, z_maxshear ,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('Max shear strain','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
% alpha(c,.5)



 