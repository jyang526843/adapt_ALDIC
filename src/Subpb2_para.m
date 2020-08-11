function [Uhat] = Subpb2_para(coordinatesFEM, elementsFEM, M, N, neumann, beta, mu, U, F, W, v, alpha, GaussPtOrder, clusterNo)


DIM = 2; NodesNumPerEle = 4; FEMSize = DIM*(size(coordinatesFEM,1));

% ====== Set boundary conditions ======
% No Dirichlet boundary conditions
dirichlet = []; % Set boundary condition by default that dirichlet = [];
% Set Neumann boundary conditions using input "F" gradient tensors.
winstepsize = abs(coordinatesFEM(1,1) - coordinatesFEM(2,1));

% ====== Initialize FE-system ======
% A = sparse(2*size(coordinatesFEM,1),2*size(coordinatesFEM,1));
% b = sparse(2*size(coordinatesFEM,1),1);

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
 

% if clusterNo== 0 || clusterNo== 1
%     hbar = waitbar(0, ['FE-method to solve Subproblem 2.']);
% else
%     hbar=parfor_progressbar(size(elementsFEM,1),['FE-method to solve Subproblem 2.']);
% end

% ============= Each element, assemble stiffness matrix ============
parfor j = 1:size(elementsFEM,1)
      
    if clusterNo== 0 || clusterNo== 1
        waitbar(j/size(elementsFEM,1));
    else
        hbar.iterate(1);
    end
    
    tempA = zeros(DIM*4,DIM*4); tempb = tempA(:,1); 
    
    % ------ Find four corner points ------
    point1x = coordinatesFEM(elementsFEM(j,1),1); point1y = coordinatesFEM(elementsFEM(j,1),2);
    point2x = coordinatesFEM(elementsFEM(j,2),1); point2y = coordinatesFEM(elementsFEM(j,2),2);
    point3x = coordinatesFEM(elementsFEM(j,3),1); point3y = coordinatesFEM(elementsFEM(j,3),2);
    point4x = coordinatesFEM(elementsFEM(j,4),1); point4y = coordinatesFEM(elementsFEM(j,4),2);
     
    % ------ Find the element nodal indices ------
    tempIndexU = [2*elementsFEM(j,1)-1 2*elementsFEM(j,1) 2*elementsFEM(j,2)-1 2*elementsFEM(j,2)...
            2*elementsFEM(j,3)-1 2*elementsFEM(j,3) 2*elementsFEM(j,4)-1 2*elementsFEM(j,4)];
    % ------ Find deformation gradient tensors indices ------
    tempIndexF = [4*elementsFEM(j,1)-3 4*elementsFEM(j,1)-1 4*elementsFEM(j,1)-2 4*elementsFEM(j,1); ...
             4*elementsFEM(j,1)-3 4*elementsFEM(j,1)-1 4*elementsFEM(j,1)-2 4*elementsFEM(j,1); ...
             4*elementsFEM(j,2)-3 4*elementsFEM(j,2)-1 4*elementsFEM(j,2)-2 4*elementsFEM(j,2); ...
             4*elementsFEM(j,2)-3 4*elementsFEM(j,2)-1 4*elementsFEM(j,2)-2 4*elementsFEM(j,2); ...
             4*elementsFEM(j,3)-3 4*elementsFEM(j,3)-1 4*elementsFEM(j,3)-2 4*elementsFEM(j,3); ...
             4*elementsFEM(j,3)-3 4*elementsFEM(j,3)-1 4*elementsFEM(j,3)-2 4*elementsFEM(j,3); ...
             4*elementsFEM(j,4)-3 4*elementsFEM(j,4)-1 4*elementsFEM(j,4)-2 4*elementsFEM(j,4);
             4*elementsFEM(j,4)-3 4*elementsFEM(j,4)-1 4*elementsFEM(j,4)-2 4*elementsFEM(j,4)];
    
    % ------ Compute DivFMinusW ------
    FMinusW = F(tempIndexF) - W(tempIndexF);
    % % DivFMinusW1 = 0.5*(  -FMinusW(4*1-3)+FMinusW(4*2-3)+FMinusW(4*3-3)-FMinusW(4*4-3) ...
    % %     -FMinusW(4*1-1)-FMinusW(4*2-1)+FMinusW(4*3-1)+FMinusW(4*4-1)  );
    % % DivFMinusW2 = 0.5*(  -FMinusW(4*1-2)+FMinusW(4*2-2)+FMinusW(4*3-2)-FMinusW(4*4-2) ...
    % %     -FMinusW(4*1)-FMinusW(4*2)+FMinusW(4*3)+FMinusW(4*4)  );
    % % DivFMinusW = [DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2];
    
    % DivFMinusW = DivFMinusWGlobal(temp);           
    
    
    % ------ Compute UMinusV ------
    UMinusV = U(tempIndexU) - v(tempIndexU); 
    
    
    % ------ Nine Gauss integral points ------
    ksietaList = [-sqrt(3/5), 0, sqrt(3/5)];
    ksietaWeightList = [5/9, 8/9, 5/9];
    
    [ksiAll, etaAll] = ndgrid(ksietaList, ksietaList);
    [weightksiAll, weightetaAll] = ndgrid(ksietaWeightList, ksietaWeightList);
    
    ksiAll = ksiAll(:); etaAll = etaAll(:); weightksiAll = weightksiAll(:); weightetaAll = weightetaAll(:);
    
    
    %for tempi = 1:3
    %    for tempj = 1:3
    for tempjj = 1:length(ksiAll)
        
        % ---------------------------
            % ksi = ksietaList(tempi);
            % eta = ksietaList(tempj);
            % weightksi = ksietaWeightList(tempi);
            % weighteta = ksietaWeightList(tempj);
            ksi = ksiAll(tempjj); eta = etaAll(tempjj); weightksi = weightksiAll(tempjj); weighteta = weightetaAll(tempjj);
            
            
            % ------ Calculate N ------
            N1 = (1-ksi)*(1-eta)*0.25;
            N2 = (1+ksi)*(1-eta)*0.25;
            N3 = (1+ksi)*(1+eta)*0.25;
            N4 = (1-ksi)*(1+eta)*0.25;
            N = [N1 0 N2 0 N3 0 N4 0;
                0 N1 0 N2 0 N3 0 N4];
            NDiag = diag([N1,N1,N2,N2,N3,N3,N4,N4]);
            
            % ------ Calculate Jacobian matrix ------
            J11 = funDN1x(ksi,eta)*point1x + funDN2x(ksi,eta)*point2x + ...
                  funDN3x(ksi,eta)*point3x + funDN4x(ksi,eta)*point4x;
            J12 = funDN1x(ksi,eta)*point1y + funDN2x(ksi,eta)*point2y + ...
                  funDN3x(ksi,eta)*point3y + funDN4x(ksi,eta)*point4y;
            J21 = funDN1y(ksi,eta)*point1x + funDN2y(ksi,eta)*point2x + ...
                  funDN3y(ksi,eta)*point3x + funDN4y(ksi,eta)*point4x;
            J22 = funDN1y(ksi,eta)*point1y + funDN2y(ksi,eta)*point2y + ...
                  funDN3y(ksi,eta)*point3y + funDN4y(ksi,eta)*point4y;
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
            %A(tempIndexU,tempIndexU) = A(tempIndexU,tempIndexU) + Jacobian*weightksi*weighteta * ( beta*(DN')*(DN) + mu*(N')*N + alpha*mu*(DN')*(DN) );
             tempA = tempA + Jacobian*weightksi*weighteta * ( beta*(DN')*(DN) + mu*(N')*N + alpha*mu*(DN')*(DN) );
            
            % ------ Construct b vector ------
            % b(temp) = b(temp) + Jacobian*weightksi*weighteta* ( -beta*NDiag*(DivFMinusW) + mu*NDiag*(UMinusV) );
            % b(tempIndexU) = b(tempIndexU) + Jacobian*weightksi*weighteta* ( beta*diag((DN')*(FMinusW')) + mu*(N')*N*(UMinusV) );
            tempb = tempb + Jacobian*weightksi*weighteta* ( beta*diag((DN')*(FMinusW')) + mu*(N')*N*(UMinusV) );
        
    %   end
    % end
    end
    
    [IndexAXX,IndexAYY] = ndgrid(tempIndexU,tempIndexU);
    INDEXAIpar{j}=IndexAXX(:); INDEXAJpar{j}=IndexAYY(:); INDEXAVALpar{j}=tempA(:);
    INDEXBIpar{j}=tempIndexU(:); INDEXBVALpar{j}=tempb(:);   
    % --- To store A_ele for each element ---
    %A_ele(j,1:64)=reshape(A(temp,temp),1,64);
    
end

% close(hbar);


% ====== Initialize A matrix and b vector ======
A = sparse(FEMSize , FEMSize);   b=sparse(FEMSize ,1);
for eleInd = 1:size(elementsFEM,1)
    A = A + sparse(INDEXAIpar{eleInd}, INDEXAJpar{eleInd}, INDEXAVALpar{eleInd}, FEMSize  , FEMSize ) ;
    b = b + sparse(INDEXBIpar{eleInd}, ones(length(INDEXBIpar{eleInd}),1),INDEXBVALpar{eleInd}, FEMSize  , 1) ;
end
    
    

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
