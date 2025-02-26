% ==============================================
% function Plotdisp_show
% ==============================================
function Plotdisp_show(U0,coordinatesFEM,elementsFEM,varargin)

switch nargin
    case 4
        edgeColorOrNot = varargin{1};
    otherwise
        edgeColorOrNot = 'EdgeColor';
end

warning off;

figure; show([],elementsFEM(:,1:4),coordinatesFEM,U0(1:2:end),edgeColorOrNot); 
title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
view(2); set(gca,'fontsize',18); axis tight; axis equal; colorbar;
% view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w'); colormap jet;
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

colormap jet; 
box on



figure; show([],elementsFEM(:,1:4),coordinatesFEM,U0(2:2:end),edgeColorOrNot); 
title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
view(2); set(gca,'fontsize',18); axis tight; axis equal; colorbar;

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w'); colormap jet;
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

colormap jet;  
box on
