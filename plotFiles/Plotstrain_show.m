% ==============================================
% function Plotstrain_show
% ==============================================
 
function Plotstrain_show(FSubpb2,coordinatesFEM,elementsFEM,varargin)

switch nargin
    case 4
        edgeColorOrNot = varargin{1};
    otherwise
        edgeColorOrNot = 'EdgeColor';
end

figure; show([],elementsFEM(:,1:4),coordinatesFEM,FSubpb2(1:4:end),edgeColorOrNot); 
title('Strain $e_{11}$','fontweight','normal','Interpreter','latex'); set(gca,'fontsize',18); 
view(2); axis tight; axis equal; 
colorbar;  box on;

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
colormap jet;


figure; show([],elementsFEM(:,1:4),coordinatesFEM,FSubpb2(4:4:end),edgeColorOrNot); 
title('Strain $e_{22}$','fontweight','normal','Interpreter','latex'); set(gca,'fontsize',18); 
view(2); axis tight; axis equal; 
colorbar;   box on;

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
colormap jet;


figure; show([],elementsFEM(:,1:4),coordinatesFEM,0.5*(FSubpb2(2:4:end)+FSubpb2(3:4:end)),edgeColorOrNot); 
title('Strain $e_{12}$','fontweight','normal','Interpreter','latex'); set(gca,'fontsize',18); 
view(2); axis tight; axis equal;
colorbar;  box on;

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
colormap jet;



F_maxshear = sqrt((0.5*(FSubpb2(1:4:end)-FSubpb2(4:4:end))).^2 + (0.5*FSubpb2(2:4:end)+0.5*FSubpb2(3:4:end)).^2);
figure; show([],elementsFEM(:,1:4),coordinatesFEM,F_maxshear,edgeColorOrNot); 
title('Max shear','fontweight','normal','Interpreter','latex'); set(gca,'fontsize',18); 
view(2); axis tight; axis equal;
colorbar;  box on;

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
colormap jet;

 


