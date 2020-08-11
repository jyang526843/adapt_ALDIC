
function Plotstrain(U,Rad,FSubpb3,x,y,sizeOfImg,CurrentImg,OrigDICImgTransparency)

warning off;

M = size(x,1); N = size(x,2);
u_x = FSubpb3(1:4:end); v_x = FSubpb3(2:4:end);
u_y = FSubpb3(3:4:end); v_y = FSubpb3(4:4:end);
 
u_x = reshape(u_x,M,N); v_x = reshape(v_x,M,N);
u_y = reshape(u_y,M,N); v_y = reshape(v_y,M,N);

u = U(1:2:end); v = U(2:2:end);
u0 = reshape(u,M+2*Rad,N+2*Rad); v0 = reshape(v,M+2*Rad,N+2*Rad);
u = u0(1+Rad:end-Rad,1+Rad:end-Rad); v = v0(1+Rad:end-Rad,1+Rad:end-Rad);

% imagesc([x(1,1) x(end,1)], [y(1,1) y(1,end)], flipud(g)); hold on;
if M < 9,x2 = x(:,1)'; else x2 = interp(x(:,1)',4); end
if N < 9,y2 = y(1,:); else y2 = interp(y(1,:),4); end



z_u = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u,M*N,1),x2,y2);

z_v = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v,M*N,1),x2,y2);


z_exx = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u_x,M*N,1),x2,y2);
z_exy = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(0.5*(v_x+u_y),M*N,1),x2,y2);
z_eyy = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v_y,M*N,1),x2,y2);
z_maxshear = sqrt((0.5*(z_exx-z_eyy)).^2 + z_exy.^2);


[x2,y2]=ndgrid(x2,y2);x2=x2'; y2=y2';


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
fig1=figure; %1
ax1=axes; h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
%h1=imshow( (imread(CurrentImg)),'InitialMagnification','fit');
axis on; axis equal; axis tight; box on;
set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');

hold on; ax2=axes; h2=surf(x2+z_u,sizeOfImg(2)+1-(y2-z_v),z_exx,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis([-0.01,0.08]); %set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight;
colormap(coolwarm(32));  colormap jet;
%%Link them together
linkaxes([ax1,ax2]);  
%%Hide the top axes
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'gray'); %colormap(ax2,coolwarmmap); 

if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
 
ax1.Visible = 'on';
ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.01 .11 .03 .815]); 
cb2.TickLabelInterpreter = 'latex';
set(gcf,'color','w');
xlabel( '$x$ (pixels)','Interpreter','latex'); 
ylabel('$y$ (pixels)','Interpreter','latex');
title('Strain $e_{xx}$','FontWeight','Normal','Interpreter','latex');
  




fig1=figure; %1
ax1=axes; h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
%h1=imshow( (imread(CurrentImg)),'InitialMagnification','fit');
axis on; axis equal; axis tight; box on;
set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');

hold on; ax2=axes; h2=surf(x2+z_u,sizeOfImg(2)+1-(y2-z_v),z_exy,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on; caxis([-0.03,0.03]); %set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight;
colormap(coolwarm(32));  colormap jet;
%%Link them together
linkaxes([ax1,ax2]);  
%%Hide the top axes
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'gray'); % colormap(ax2,coolwarmmap); 

if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
 
ax1.Visible = 'on';
ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.01 .11 .03 .815]); 
cb2.TickLabelInterpreter = 'latex';
set(gcf,'color','w');
xlabel( '$x$ (pixels)','Interpreter','latex'); 
ylabel('$y$ (pixels)','Interpreter','latex');
title('Strain $e_{xy}$','FontWeight','Normal','Interpreter','latex');



fig1=figure; %1
ax1=axes; h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
%h1=imshow( (imread(CurrentImg)),'InitialMagnification','fit');
axis on; axis equal; axis tight; box on;
set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');

hold on; ax2=axes; h2=surf(x2+z_u,sizeOfImg(2)+1-(y2-z_v),z_eyy,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis([-0.08,0.01]); %set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight;
colormap(coolwarm(32));  colormap jet;
%%Link them together
linkaxes([ax1,ax2]);  
%%Hide the top axes
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'gray'); % colormap(ax2,coolwarmmap); 

if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
 
ax1.Visible = 'on';
ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.01 .11 .03 .815]); 
cb2.TickLabelInterpreter = 'latex';
set(gcf,'color','w');
xlabel( '$x$ (pixels)','Interpreter','latex'); 
ylabel('$y$ (pixels)','Interpreter','latex');
title('Strain $e_{yy}$','FontWeight','Normal','Interpreter','latex');



  
fig1=figure; %1
ax1=axes; h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
%h1=imshow( (imread(CurrentImg)),'InitialMagnification','fit');
axis on; axis equal; axis tight; box on;
set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');

hold on; ax2=axes; h2=surf(x2+z_u,sizeOfImg(2)+1-(y2-z_v),z_maxshear,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis([0,0.07]); %set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight;
colormap(coolwarm(32));  colormap jet;
%%Link them together
linkaxes([ax1,ax2]);  
%%Hide the top axes
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'gray'); % colormap(ax2,coolwarmmap);  colormap( coolwarmmap);

if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
 
ax1.Visible = 'on';
ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.01 .11 .03 .815]); 
cb2.TickLabelInterpreter = 'latex';
set(gcf,'color','w');
xlabel( '$x$ (pixels)','Interpreter','latex'); 
ylabel('$y$ (pixels)','Interpreter','latex');
title('Max shear strain','FontWeight','Normal','Interpreter','latex');




end



% function Plotstrain(FSubpb3,x,y,f,g)
%   
% M = size(x,1); N = size(x,2);
% u_x = FSubpb3(1:4:end); v_x = FSubpb3(2:4:end);
% u_y = FSubpb3(3:4:end); v_y = FSubpb3(4:4:end);
%  
% u_x = reshape(u_x,M,N); v_x = reshape(v_x,M,N);
% u_y = reshape(u_y,M,N); v_y = reshape(v_y,M,N);
% 
% %%%%%%%%%%% Plot of exx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  
%  figure; 
% fig = surf(uint8(rot90(f,1)),'EdgeColor','none','LineStyle','none');
% axis tight; axis equal; axis off; view(2); set(gcf,'color','w');
% xl = xlim ; yl = ylim;
% colormap gray
% set(gca,'fontSize',18); 
% b = colorbar;  b.TickLabelInterpreter = 'latex'; delete(b)
% 
% drawnow
% ax = gca;
% ax.Units = 'pixels';
% pos = ax.Position;
%  
% rect = [ 46, 0, pos(3)-92, pos(4)];
% F = getframe(gca,rect);
% im0 = F.cdata;
% % im0 = print('-RGBImage','-painters','-r600'); %print('-RGBImage','-r600');
% 
% clf;
% if M < 9 
%     x2 = x(:,1)';
% else  
%     x2 = interp(x(:,1)',4);
% end
% if N < 9
%     y2 = y(1,:);
% else
% 	y2 = interp(y(1,:),4);
% end
% z = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u_x,M*N,1),x2,y2);
% % [c,h] = contourf(x2,y2,z,20,'EdgeColor','none','LineStyle','none');
% surf(x2,y2,z,'EdgeColor','none','LineStyle','none')
% axis equal; axis tight; axis off; set(gca,'fontSize',18); view(2);
% set(gcf,'color','none'); xlim(xl); ylim(yl);  
% colormap jet
% 
% b = colorbar; caxis([-0.35,0.4]); b.TickLabelInterpreter = 'latex'; caxisGet = caxis; delete(b);
% 
% if x(M,N) < 200
%     set(gca,'XTick',[]);
% end
% if y(M,N) < 200
%     set(gca,'YTick',[]);
% end
%  
% drawnow
% ax = gca;
% ax.Units = 'pixels';
% pos = ax.Position;
%  
% rect = [46, 0, pos(3)-92, pos(4)];
% F = getframe(gca,rect);
% im1 = F.cdata;
% 
% clf
% h0 = imshow(im0); 
% 
% hold on
% h1 = imshow(im1);
% set(h1,'AlphaData',0.7);
% 
% set(gca,'fontSize',18);
% 
% title('Strain $e_{xx}$','FontWeight','Normal','Interpreter','latex');
% xlim([0,pos(3)-92]); ylim([0,pos(4)]);
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% % xlabh = get(gca,'xlabel'); set(xlabh,'Position',get(xlabh,'Position') - [0 15 0]);
% % ylabh = get(gca,'ylabel'); set(ylabh,'Position',get(ylabh,'Position') + [25 0 0]);
% 
% set(gcf,'color','w');
% a = gca; a.TickLabelInterpreter = 'latex';
% colormap jet; 
% b = colorbar; b.TickLabelInterpreter = 'latex'; caxis(caxisGet);  set(gca,'fontSize',18);
% box on
% 
% axis on
% xlength = pos(3)-92; 
% xRatio = size(g,1)/xlength;
% XTickOrig = [0, 100, 200, 300, 400, 500];
% XTickNow = XTickOrig/xRatio;
% ax = gca;
% ax.XTick =  XTickNow;
% xticklabels = get(ax,'XTickLabel');
% xticklabels{1} = '0';
% xticklabels{2} = '100';
% xticklabels{3} = '200';
% xticklabels{4} = '300';
% xticklabels{5} = '400';
% xticklabels{6} = '500';
% set(ax,'XTickLabel',xticklabels);
% 
% ylength = pos(4);
% yRatio = size(g,1)/ylength;
% YTickOrig = [0, 100, 200, 300, 400, 500];
% YTickNow = YTickOrig/yRatio;
% ax = gca;
% ax.YTick =  YTickNow;
% yticklabels = get(ax,'YTickLabel');
% yticklabels{1} = '0';
% yticklabels{2} = '100';
% yticklabels{3} = '200';
% yticklabels{4} = '300';
% yticklabels{5} = '400';
% yticklabels{6} = '500';
% set(ax,'YTickLabel',yticklabels);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%% Plot of eyy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  
% figure; 
% fig = surf(uint8(rot90(f,1)),'EdgeColor','none','LineStyle','none');
% axis tight; axis equal; axis off; view(2); set(gcf,'color','w');
% xl = xlim ; yl = ylim;
% colormap gray
% set(gca,'fontSize',18); 
% b = colorbar; b.TickLabelInterpreter = 'latex'; delete(b)
% 
% drawnow
% ax = gca;
% ax.Units = 'pixels';
% pos = ax.Position;
%  
% rect = [ 46, 0, pos(3)-92, pos(4)];
% F = getframe(gca,rect);
% im0 = F.cdata;
% % im0 = print('-RGBImage','-painters','-r600'); %print('-RGBImage','-r600');
% 
% clf;
% if M < 9 
%     x2 = x(:,1)';
% else  
%     x2 = interp(x(:,1)',4);
% end
% if N < 9
%     y2 = y(1,:);
% else
% 	y2 = interp(y(1,:),4);
% end
% z = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v_y,M*N,1),x2,y2);
% % [c,h] = contourf(x2,y2,z,20,'EdgeColor','none','LineStyle','none');
% surf(x2,y2,z,'EdgeColor','none','LineStyle','none')
% axis equal; axis tight; axis off; set(gca,'fontSize',18); view(2);
% set(gcf,'color','none'); xlim(xl); ylim(yl);  
% colormap jet
% 
% b = colorbar; caxis([-0.5,0.4]); b.TickLabelInterpreter = 'latex'; caxisGet = caxis; delete(b);
% 
% if x(M,N) < 200
%     set(gca,'XTick',[]);
% end
% if y(M,N) < 200
%     set(gca,'YTick',[]);
% end
%  
% drawnow
% ax = gca;
% ax.Units = 'pixels';
% pos = ax.Position;
%  
% rect = [46, 0, pos(3)-92, pos(4)];
% F = getframe(gca,rect);
% im1 = F.cdata;
% 
% clf
% h0 = imshow(im0); 
% 
% hold on
% h1 = imshow(im1);
% set(h1,'AlphaData',0.7);
% 
% set(gca,'fontSize',18);
% 
% title('Strain $e_{yy}$','FontWeight','Normal','Interpreter','latex');
% xlim([0,pos(3)-92]); ylim([0,pos(4)]);
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% % xlabh = get(gca,'xlabel'); set(xlabh,'Position',get(xlabh,'Position') - [0 15 0]);
% % ylabh = get(gca,'ylabel'); set(ylabh,'Position',get(ylabh,'Position') + [25 0 0]);
% 
% set(gcf,'color','w');
% a = gca; a.TickLabelInterpreter = 'latex';
% colormap jet; 
% b = colorbar; b.TickLabelInterpreter = 'latex'; caxis(caxisGet);  set(gca,'fontSize',18);
% box on
% 
% axis on
% xlength = pos(3)-92; 
% xRatio = size(g,1)/xlength;
% XTickOrig = [0, 100, 200, 300, 400, 500];
% XTickNow = XTickOrig/xRatio;
% ax = gca;
% ax.XTick =  XTickNow;
% xticklabels = get(ax,'XTickLabel');
% xticklabels{1} = '0';
% xticklabels{2} = '100';
% xticklabels{3} = '200';
% xticklabels{4} = '300';
% xticklabels{5} = '400';
% xticklabels{6} = '500';
% set(ax,'XTickLabel',xticklabels);
% 
% ylength = pos(4);
% yRatio = size(g,1)/ylength;
% YTickOrig = [0, 100, 200, 300, 400, 500];
% YTickNow = YTickOrig/yRatio;
% ax = gca;
% ax.YTick =  YTickNow;
% yticklabels = get(ax,'YTickLabel');
% yticklabels{1} = '0';
% yticklabels{2} = '100';
% yticklabels{3} = '200';
% yticklabels{4} = '300';
% yticklabels{5} = '400';
% yticklabels{6} = '500';
% set(ax,'YTickLabel',yticklabels);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% %%%%%%%%%%% Plot of exy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure; 
% fig = surf(uint8(rot90(f,1)),'EdgeColor','none','LineStyle','none');
% axis tight; axis equal; axis off; view(2); set(gcf,'color','w');
% xl = xlim ; yl = ylim;
% colormap gray
% set(gca,'fontSize',18); 
% b = colorbar; b.TickLabelInterpreter = 'latex'; delete(b)
% 
% drawnow
% ax = gca;
% ax.Units = 'pixels';
% pos = ax.Position;
%  
% rect = [ 46, 0, pos(3)-92, pos(4)];
% F = getframe(gca,rect);
% im0 = F.cdata;
% % im0 = print('-RGBImage','-painters','-r600'); %print('-RGBImage','-r600');
% 
% clf;
% if M < 9 
%     x2 = x(:,1)';
% else  
%     x2 = interp(x(:,1)',4);
% end
% if N < 9
%     y2 = y(1,:);
% else
% 	y2 = interp(y(1,:),4);
% end
% z = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(0.5*(u_x+v_y),M*N,1),x2,y2);
% % [c,h] = contourf(x2,y2,z,20,'EdgeColor','none','LineStyle','none');
% surf(x2,y2,z,'EdgeColor','none','LineStyle','none')
% axis equal; axis tight; axis off; set(gca,'fontSize',18); view(2);
% set(gcf,'color','none'); xlim(xl); ylim(yl);  
% colormap jet
% 
% b = colorbar; caxis([-0.4,0.4]); b.TickLabelInterpreter = 'latex'; caxisGet = caxis; delete(b);
% 
% if x(M,N) < 200
%     set(gca,'XTick',[]);
% end
% if y(M,N) < 200
%     set(gca,'YTick',[]);
% end
%  
% drawnow
% ax = gca;
% ax.Units = 'pixels';
% pos = ax.Position;
%  
% rect = [46, 0, pos(3)-92, pos(4)];
% F = getframe(gca,rect);
% im1 = F.cdata;
% 
% clf
% h0 = imshow(im0); 
% 
% hold on
% h1 = imshow(im1);
% set(h1,'AlphaData',0.7);
% 
% set(gca,'fontSize',18);
% 
% title('Strain $e_{xy}$','FontWeight','Normal','Interpreter','latex');
% xlim([0,pos(3)-92]); ylim([0,pos(4)]);
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% % xlabh = get(gca,'xlabel'); set(xlabh,'Position',get(xlabh,'Position') - [0 15 0]);
% % ylabh = get(gca,'ylabel'); set(ylabh,'Position',get(ylabh,'Position') + [25 0 0]);
% 
% set(gcf,'color','w');
% a = gca; a.TickLabelInterpreter = 'latex';
% colormap jet; 
% b = colorbar; b.TickLabelInterpreter = 'latex'; caxis(caxisGet);  set(gca,'fontSize',18);
% box on
% 
% axis on
% xlength = pos(3)-92; 
% xRatio = size(g,1)/xlength;
% XTickOrig = [0, 100, 200, 300, 400, 500];
% XTickNow = XTickOrig/xRatio;
% ax = gca;
% ax.XTick =  XTickNow;
% xticklabels = get(ax,'XTickLabel');
% xticklabels{1} = '0';
% xticklabels{2} = '100';
% xticklabels{3} = '200';
% xticklabels{4} = '300';
% xticklabels{5} = '400';
% xticklabels{6} = '500';
% set(ax,'XTickLabel',xticklabels);
% 
% ylength = pos(4);
% yRatio = size(g,1)/ylength;
% YTickOrig = [0, 100, 200, 300, 400, 500];
% YTickNow = YTickOrig/yRatio;
% ax = gca;
% ax.YTick =  YTickNow;
% yticklabels = get(ax,'YTickLabel');
% yticklabels{1} = '0';
% yticklabels{2} = '100';
% yticklabels{3} = '200';
% yticklabels{4} = '300';
% yticklabels{5} = '400';
% yticklabels{6} = '500';
% set(ax,'YTickLabel',yticklabels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 