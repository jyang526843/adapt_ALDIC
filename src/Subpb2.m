function [Uhat] = Subpb2(coordinatesFEM, elementsFEM, M, N, neumann, beta, mu, U, F, W, v, alpha, GaussPtOrder, clusterNo)

DIM = 2; NodesPerEle = 4; % Using cubic elements
FEMSize = DIM*size(coordinatesFEM,1); % FE-system size

% ====== Set boundary conditions ======
% No Dirichlet boundary conditions
dirichlet = []; % Set boundary condition by default that dirichlet = [];
% Set Neumann boundary conditions using input "F" gradient tensors.
winstepsize = abs(coordinatesFEM(1,1) - coordinatesFEM(2,1));

% ====== Initialize FE-system ======
%A = sparse(2*size(coordinatesFEM,1),2*size(coordinatesFEM,1));
%b = sparse(2*size(coordinatesFEM,1),1);
INDEXAI = []; INDEXAJ = []; INDEXAVAL = []; 
INDEXBI = []; INDEXBVAL = [];

% ====== Compute Div(F-W) ======
% FMinusW = F - W; 
% FMinusW11 = reshape(FMinusW(1:4:end),M,N); FMinusW21 = reshape(FMinusW(2:4:end),M,N);
% FMinusW12 = reshape(FMinusW(3:4:end),M,N); FMinusW22 = reshape(FMinusW(4:4:end),M,N);

% imgradientMatrix = [-1/2 0 1/2]'; % Central finite difference operator
% DFMinusW11Dxtemp = imfilter(FMinusW11,imgradientMatrix);
% DFMinusW12Dytemp = imfilter(FMinusW12,imgradientMatrix');
% DFMinusW11Dx     = DFMinusW11Dxtemp(1:end, 1:end);
% DFMinusW12Dy     = DFMinusW12Dytemp(1:end, 1:end);
% 
% DFMinusW21Dxtemp = imfilter(FMinusW21,imgradientMatrix);
% DFMinusW22Dytemp = imfilter(FMinusW22,imgradientMatrix');
% DFMinusW21Dx     = DFMinusW21Dxtemp(1:end, 1:end);
% DFMinusW22Dy     = DFMinusW22Dytemp(1:end, 1:end);
% 
% DivFMinusW1temp = DFMinusW11Dx + DFMinusW12Dy;
% DivFMinusW2temp = DFMinusW21Dx + DFMinusW22Dy;
% 
% DivFMinusWGlobal  = zeros(M*N*2,1);
% DivFMinusWGlobal(1:2:end) = reshape(DivFMinusW1temp,M*N,1);
% DivFMinusWGlobal(2:2:end) = reshape(DivFMinusW2temp,M*N,1);
% 
% DFMinusW11DxStartx = 1; DFMinusW11DyStarty = 1;
 
hbar = waitbar(0,['Global step']);
% ============= Each element, assemble stiffness matrix ============
for j = 1:size(elementsFEM,1)
     
    waitbar(j/size(elementsFEM,1));
    
    tempA = zeros(DIM*NodesPerEle,DIM*NodesPerEle); tempb = tempA(:,1);
        
    % ------ Find four corner pts ------
    pt1xyz = coordinatesFEM(elementsFEM(j,1),:); pt2xyz = coordinatesFEM(elementsFEM(j,2),:);
    pt3xyz = coordinatesFEM(elementsFEM(j,3),:); pt4xyz = coordinatesFEM(elementsFEM(j,4),:);
    
    pt1x = pt1xyz(1); pt1y = pt1xyz(2);  pt2x = pt2xyz(1); pt2y = pt2xyz(2);
    pt3x = pt3xyz(1); pt3y = pt3xyz(2);  pt4x = pt4xyz(1); pt4y = pt4xyz(2);
         
    % ------ Find the element nodal indices ------
    %temp = [2*elementsFEM(j,1)-1 2*elementsFEM(j,1) 2*elementsFEM(j,2)-1 2*elementsFEM(j,2)...
    %        2*elementsFEM(j,3)-1 2*elementsFEM(j,3) 2*elementsFEM(j,4)-1 2*elementsFEM(j,4)];
    tp = ones(1,DIM);
    tempIndexU = 2*elementsFEM(j,[tp,2*tp,3*tp,4*tp]);
    tempIndexU(1:2:end) = tempIndexU(1:2:end)-1;
    
    % ------ Find deformation gradient tensors indices ------
    tempF = [4*elementsFEM(j,1)-3 4*elementsFEM(j,1)-1 4*elementsFEM(j,1)-2 4*elementsFEM(j,1); ...
             4*elementsFEM(j,1)-3 4*elementsFEM(j,1)-1 4*elementsFEM(j,1)-2 4*elementsFEM(j,1); ...
             4*elementsFEM(j,2)-3 4*elementsFEM(j,2)-1 4*elementsFEM(j,2)-2 4*elementsFEM(j,2); ...
             4*elementsFEM(j,2)-3 4*elementsFEM(j,2)-1 4*elementsFEM(j,2)-2 4*elementsFEM(j,2); ...
             4*elementsFEM(j,3)-3 4*elementsFEM(j,3)-1 4*elementsFEM(j,3)-2 4*elementsFEM(j,3); ...
             4*elementsFEM(j,3)-3 4*elementsFEM(j,3)-1 4*elementsFEM(j,3)-2 4*elementsFEM(j,3); ...
             4*elementsFEM(j,4)-3 4*elementsFEM(j,4)-1 4*elementsFEM(j,4)-2 4*elementsFEM(j,4);
             4*elementsFEM(j,4)-3 4*elementsFEM(j,4)-1 4*elementsFEM(j,4)-2 4*elementsFEM(j,4)];
    
    % ------ Compute DivFMinusW ------
    FMinusW = F(tempF) - W(tempF);
    % % DivFMinusW1 = 0.5*(  -FMinusW(4*1-3)+FMinusW(4*2-3)+FMinusW(4*3-3)-FMinusW(4*4-3) ...
    % %     -FMinusW(4*1-1)-FMinusW(4*2-1)+FMinusW(4*3-1)+FMinusW(4*4-1)  );
    % % DivFMinusW2 = 0.5*(  -FMinusW(4*1-2)+FMinusW(4*2-2)+FMinusW(4*3-2)-FMinusW(4*4-2) ...
    % %     -FMinusW(4*1)-FMinusW(4*2)+FMinusW(4*3)+FMinusW(4*4)  );
    % % DivFMinusW = [DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2];
    
    % DivFMinusW = DivFMinusWGlobal(tempIndexU);           
    
    % ------ Compute UMinusV ------
    UMinusV = U(tempIndexU) - v(tempIndexU); 
    
    % ------ Nine Gauss integral pts ------
    ksietaList = [-sqrt(3/5), 0, sqrt(3/5)];
    ksietaWeightList = [5/9, 8/9, 5/9];
    [ksiAll,etaAll] = ndgrid(ksietaList,ksietaList);
    [ksiWeightAll,etaWeightAll] = ndgrid(ksietaWeightList,ksietaWeightList);
    ksiAll = ksiAll(:); etaAll = etaAll(:);
    ksiWeightAll = ksiWeightAll(:); etaWeightAll = etaWeightAll(:);
    
    %for tempi = 1:3
    %    for tempj = 1:3
    for tempjj = 1:length(ksiAll)
        
            ksi = ksiAll(tempjj);
            eta = etaAll(tempjj);
            weightksi = ksiWeightAll(tempjj);
            weighteta = etaWeightAll(tempjj);
            
            % ------ Calculate N ------
            N1 = (1-ksi)*(1-eta)*0.25;
            N2 = (1+ksi)*(1-eta)*0.25;
            N3 = (1+ksi)*(1+eta)*0.25;
            N4 = (1-ksi)*(1+eta)*0.25;
            N = [N1 0 N2 0 N3 0 N4 0;
                0 N1 0 N2 0 N3 0 N4];
            NDiag = diag([N1,N1,N2,N2,N3,N3,N4,N4]);
            
            % ------ Calculate Jacobian matrix ------
            J11 = funDN1x(ksi,eta)*pt1x + funDN2x(ksi,eta)*pt2x + ...
                  funDN3x(ksi,eta)*pt3x + funDN4x(ksi,eta)*pt4x;
            J12 = funDN1x(ksi,eta)*pt1y + funDN2x(ksi,eta)*pt2y + ...
                  funDN3x(ksi,eta)*pt3y + funDN4x(ksi,eta)*pt4y;
            J21 = funDN1y(ksi,eta)*pt1x + funDN2y(ksi,eta)*pt2x + ...
                  funDN3y(ksi,eta)*pt3x + funDN4y(ksi,eta)*pt4x;
            J22 = funDN1y(ksi,eta)*pt1y + funDN2y(ksi,eta)*pt2y + ...
                  funDN3y(ksi,eta)*pt3y + funDN4y(ksi,eta)*pt4y;
            J = [J11 J12; J21 J22];
            
            Jacobian = det(J);
            InvJ = 1/Jacobian*[J22 -J12; -J21 J11];
            
            % ------ DN matrix ------
            DN = [InvJ zeros(2,2); zeros(2,2) InvJ] * ...
                 [funDN1x(ksi,eta) 0 funDN2x(ksi,eta) 0 funDN3x(ksi,eta) 0 funDN4x(ksi,eta) 0;
                 funDN1y(ksi,eta) 0 funDN2y(ksi,eta) 0 funDN3y(ksi,eta) 0 funDN4y(ksi,eta) 0;
                 0 funDN1x(ksi,eta) 0 funDN2x(ksi,eta) 0 funDN3x(ksi,eta) 0 funDN4x(ksi,eta);
                 0 funDN1y(ksi,eta) 0 funDN2y(ksi,eta) 0 funDN3y(ksi,eta) 0 funDN4y(ksi,eta)];
            
            % ------ Construct A matrix ------
            % A(temp,temp) = A(temp,temp) + Jacobian*weightksi*weighteta * ( beta*(DN')*(DN) + mu*(N')*N + alpha*mu*(DN')*(DN) );
            tempA = tempA + Jacobian*weightksi*weighteta * ( beta*(DN')*(DN) + mu*(N')*N + alpha*mu*(DN')*(DN) );
            
            % ------ Construct b vector ------
            % b(temp) = b(temp) + Jacobian*weightksi*weighteta* ( -beta*NDiag*(DivFMinusW) + mu*NDiag*(UMinusV) );
            %b(temp) = b(temp) + Jacobian*weightksi*weighteta* ( beta*diag((DN')*(FMinusW')) + mu*(N')*N*(UMinusV) );
            tempb = tempb + Jacobian*weightksi*weighteta* ( beta*diag((DN')*(FMinusW')) + mu*(N')*N*(UMinusV) );
        
        %end
    end
    
    % --- To store A_ele for each element ---
    %A_ele(j,1:64)=reshape(A(temp,temp),1,64);
    [IndexAXX,IndexAYY] = ndgrid(tempIndexU,tempIndexU); 
    INDEXAI = [INDEXAI;IndexAXX(:)]; INDEXAJ = [INDEXAJ;IndexAYY(:)]; INDEXAVAL = [INDEXAVAL;tempA(:)];
    INDEXBI = [INDEXBI;tempIndexU(:)]; INDEXBVAL = [INDEXBVAL;tempb(:)];
    
end

close(hbar);


A = sparse(INDEXAI,INDEXAJ,INDEXAVAL,FEMSize,FEMSize);
b = sparse(INDEXBI,ones(length(INDEXBI),1),INDEXBVAL,FEMSize,1);


% ========= Dirichlet boundary conditions ==========
Uhat = sparse(2*size(coordinatesFEM,1),1);
Uhat(2*unique(dirichlet)) = U(2*unique(dirichlet));
Uhat(2*unique(dirichlet)-1) = U(2*unique(dirichlet)-1);
b = b - A * Uhat;

dirichlettemp = [2*dirichlet; 2*dirichlet-1];
FreeNodes = setdiff(1:2*size(coordinatesFEM,1),unique(dirichlettemp));

% ========= Neumann boundary conditions ==========
% Get boundary traction force using displacement input "U"
BCForce = -1/(winstepsize)*F; % F is the input of Subpb2.m, which is the Subproblem 1 solved gradient def tensors;

% If F solved by Subproblem 1 is too noisy at boundary, try the following code:
% [BCForce] = -1/(winstepsize*2)*funGlobal_NodalStrainAvg(coordinatesFEM,elementsFEM,U,GaussPtOrder);
for tempj = 1:size(neumann,1)
    % f1 = F11*n1 + F12*n2
    b(2*neumann(tempj,1:2)-1) = b(2*neumann(tempj,1:2)-1) + 0.5*norm(coordinatesFEM(neumann(tempj,1),:)-coordinatesFEM(neumann(tempj,2),:)) ...
        *( ( BCForce(4*neumann(tempj,1:2)-3) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)-1) * neumann(tempj,4) ) );
    % f2 = F21*n1 + F22*n2
    b(2*neumann(tempj,1:2))   = b(2*neumann(tempj,1:2)) + 0.5*norm(coordinatesFEM(neumann(tempj,1),:)-coordinatesFEM(neumann(tempj,2),:)) ...
        *( ( BCForce(4*neumann(tempj,1:2)-2) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)) * neumann(tempj,4) ) ) ;
      
    %b(2*neumann(tempj,1:2)-1) = b(2*neumann(tempj,1:2)-1) + 0.5*1 ...
    %    *( ( BCForce(4*neumann(tempj,1:2)-3) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)-1) * neumann(tempj,4) ) );
    %b(2*neumann(tempj,1:2))   = b(2*neumann(tempj,1:2)) + 0.5*1 ...
    %    *( ( BCForce(4*neumann(tempj,1:2)-2) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)) * neumann(tempj,4) ) ) ;

end

Uhat(FreeNodes) = A(FreeNodes,FreeNodes) \ b(FreeNodes);


% ========= subroutines for  FEM Q4 shape function derivatives ========

function DN1x=funDN1x(ksi,eta)
DN1x = -(1-eta)/4 ;

function DN1y=funDN1y(ksi,eta)
DN1y =  -(1-ksi)/4 ;

function DN2x=funDN2x(ksi,eta)
DN2x =  (1-eta)/4 ;

function DN2y=funDN2y(ksi,eta)
DN2y =  -(1+ksi)/4 ;

function DN3x=funDN3x(ksi,eta)
DN3x = (1+eta)/4 ;

function DN3y=funDN3y(ksi,eta)
DN3y =  (1+ksi)/4 ;

function DN4x=funDN4x(ksi,eta)
DN4x = -(1+eta)/4 ;

function DN4y=funDN4y(ksi,eta)
DN4y = (1-ksi)/4 ;
