% -----------------------------
% ALDIC with crack
% Code: yangjin@caltech.edu
% -----------------------------

%% Section 1
% ====== Clear MATLAB environment & mex set up Spline interpolation ======
close all; clear; clc; clearvars -global
fprintf('------------ Section 1 Start ------------ \n')
addpath('./func','./func_adapt','./plotFiles','./src','./imgDemo','./imgFolder','./func_adapt/refinement');
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')
% % cd("./Splines_interp/lib_matlab"); CompileLib; cd("../../");  % % mex bi-cubic spline interpolations
% % addpath("./Splines_interp/lib_matlab");
try mex -O ba_interp2.cpp; catch; end  % dbstop if error % % Old version codes.

% % ============== Global parameter initialition ===============
% thetaDorfler = 0.9; thetaDorflerList = zeros(4,1); % For mark elements needed to be refined
% DispFilterStd = 0; DispFilterSize = 0; StrainFilterStd = 0; StrainFilterSize = 0;
fprintf('------------ Section 1 Done ------------ \n \n')


%% Section 2
fprintf('------------ Section 2 Start ------------ \n')
% ====== Read images ======
[file_name,Img,winsize,winstepsize,gridxyROIRange,LoadImgMethod] = ReadImage; close all;
CrackOrNot = 0; CrackPath1 = [0,0]; CrackPath2 = [0,0]; CrackTipOrNot = 0; CrackTip = [0,0]; 
EnrHAndTipEleIndex = []; EnrTipEleIndex = []; dirichlet = []; neumann= [];
% ====== Cncomment the behind line and change the value you want ======
% gridxROIRange = [gridxROIRange1, gridxROIRange2]; gridyROIRange = [gridyROIRange1, gridyROIRange2];
% ====== Normalize images ======
[ImgNormalized,gridxyROIRange] = funNormalizeImg(Img,gridxyROIRange);  
fprintf('------------ Section 2 Done ------------ \n \n')

  
%%%%%%%%%%%%%% Deal with image sequences%%%%%%%%%%%%%%%%
if LoadImgMethod == 0 || LoadImgMethod == 1
fprintf('--- Extensions:  Additional parameters when dealing with DIC image sequence --- \n')
end
% If they are first two image frames, use above FFT searched initial guess;
% Else using the results of last image frame as the initial guess for the
% new image frame.
ResultDisp = cell(length(ImgNormalized)-1,1); 
ResultDefGrad = cell(length(ImgNormalized)-1,1);
ResultStrain = cell(length(ImgNormalized)-1,1);
if LoadImgMethod == 0 || LoadImgMethod == 1 % Image folder || Prefix of image names
    fprintf('Since we are dealing with image sequences, for each new frame,   \n')
    fprintf('do we use last frame result as initial guess or   \n')
    fprintf('Redo FFT initial guess for every new frame? \n    0: Use last frame; \n    1: Redo initial guess.  \n')
    prompt = 'Input here: '; StillFFTSearch = input(prompt);
    fprintf('\n');
    fprintf('Method to solve ALDIC Subproblem 2:    \n')
    fprintf('    1: Finite difference   \n    2: Finite element method  \n')
    prompt = 'Input here: ';
    Subpb2FDOrFEM = input(prompt);
else
    StillFFTSearch = 1;
    Subpb2FDOrFEM = 0;
end
 
fprintf('\n'); disp('--- Set up Parallel pool ---'); 
% ------ Assign parpool cluster No ------
prompt = 'How many parallel pools to open? (Put in 1 if no parallel computing): ';  ClusterNo = input(prompt);
% if ClusterNo > 1
%     delete(gcp); myCluster = parcluster('local'); delete(myCluster.Jobs);
%     parpool(ClusterNo,'SpmdEnabled',false);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for ImgSeqNum = 2:length(ImgNormalized)
     
    %% Section 3
    fprintf('\n'); fprintf('------------ Section 3 Start ------------ \n')
    % % Do first deformed image at first.
    gridxROIRange = gridxyROIRange.gridx; gridyROIRange = gridxyROIRange.gridy; 
    fNormalized = ImgNormalized{1}; gNormalized = ImgNormalized{ImgSeqNum};
    
    if ImgSeqNum == 2 || StillFFTSearch == 1
        % ====== Integer Search ======
        [SizeOfFFTSearchRegion,x0,y0,u,v,cc]= IntegerSearch(gridxROIRange,gridyROIRange,fNormalized,gNormalized,winsize,winstepsize,file_name);
        % ====== FEM mesh set up ======
        [coordinatesFEM,elementsFEM,coordinates,elements,dirichlet,neumann,x0,y0] = MeshSetUp(x0,y0,winstepsize);
        dirichlet=[]; 
        % ====== Remove outliers ======
        % Already included in the above IntegerSearch.m function
        % qDICOrNot = 0; Thr0 = 100; [u,v,cc] = funRemoveOutliers(u,v,cc,qDICOrNot,Thr0);
        % ====== Initial Value ======
        U0 = Init(u,v,cc.max,x0,y0,0); PlotuvInit; 
    else
        U0 = ResultDisp{ImgSeqNum-2}.U;
    end
    %Plotdisp_show(U0,[coordinatesFEM(:,1),size(fNormalized,2)+1-coordinatesFEM(:,2)],elementsFEM); % Plot initial values
    % ====== Get images grayscale gradients ======
    Df = funImgGradient(fNormalized,gNormalized); % % using finite difference;
    % ====== Compute f(X)-g(x+u) ======
    % PlotImgDiff(x0,y0,u,v,fNormalized,gNormalized);
    fprintf('------------ Section 3 Done ------------ \n \n')
    
    %% Section 4: FEM-based DIC 
    
    % =========== Solve: Level1 uniform mesh FEM Global DIC ============
    tol=1e-4; alphavar = 1; GaussPtOrder = 2; LevelNo = 1; tic; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ULevel1, normOfW, TimeICGN] = funGlobalICGN(coordinatesFEM,elementsFEM,Df,fNormalized,gNormalized,U0,alphavar,tol,ClusterNo);
    [FLevel1] = funGlobal_NodalStrainAvg(coordinatesFEM,elementsFEM,ULevel1,GaussPtOrder);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['--- Done! step ',num2str(LevelNo),' ---']);
    elementsFEM = [elementsFEM, zeros(size(elementsFEM,1),4)]; eleGeneration = ones(size(elementsFEM,1),1);
    LevelNo = 1; UIter = full(ULevel1); FIter = full(FLevel1); 
     
    % ============== Assign values to variables ===============
    normOfWIter = normOfW; alphaIter = alphavar;
    elementsFEMIter = elementsFEM; coordinatesFEMIter = coordinatesFEM;  
    dirichletIter = dirichlet; neumannIter = neumann; eleGenerationIter = eleGeneration;
    irregularEdgeIter = zeros(1,3);
    PassCrackOrNotIter = zeros(size(elementsFEMIter,1),1);
    try 
        if dirichletIter(1) == 0, dirichletIter = dirichletIter(2:end,1:2);
         else, dirichletIter = dirichletIter(1:end,1:2); end
    catch
    end
    % -------- Plot disp and strain field -------- 
    close all;
    Plotdisp_show(UIter,coordinatesFEMIter,elementsFEMIter(:,1:4));
    Plotstrain_show(FIter,coordinatesFEMIter,elementsFEMIter(:,1:4));
    
    figure(1); savefig(['fig_L',num2str(LevelNo),'_dispu.fig']);
    figure(2); savefig(['fig_L',num2str(LevelNo),'_dispv.fig']);
    figure(3); savefig(['fig_L',num2str(LevelNo),'_e11.fig']);
    figure(4); savefig(['fig_L',num2str(LevelNo),'_e22.fig']);
    figure(5); savefig(['fig_L',num2str(LevelNo),'_e12.fig']);
    figure(6); savefig(['fig_L',num2str(LevelNo),'_eshear.fig']);
    
    
    % ============== Store solution ===============
    solnIter = struct('UIter',UIter,'FIter',FIter,'alphaIter',alphaIter,'normOfWIter',normOfWIter, ...
        'coordinatesFEMIter',coordinatesFEMIter,'elementsFEMIter',elementsFEMIter, ...
        'eleGenerationIter',eleGenerationIter, 'irregularEdgeIter', irregularEdgeIter, ...
        'dirichletIter',dirichletIter,'neumannIter',neumannIter, ...
        'TimeICGNIter',TimeICGN,'EstimateTimeLastIter',[0],'MarkTimeLastIter',[0],'RefineTimeLastIter',[0],...
        'rhoLastIterVector',[0],'rhoLastIterVectortemp1',[0],'rhoLastIterVectortemp2',[0],'thetaDorflerLastIter',[0],...
        'EnrHAndTipEleIndex',EnrHAndTipEleIndex,'EnrTipEleIndex',EnrTipEleIndex);
    
    soln{LevelNo} = solnIter;
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while (LevelNo<6) && (winstepsize/(2^(LevelNo-1))>1) && (mod((winstepsize/(2^(LevelNo-1))),2)==0)
      
        %% Adaptive-Refine - Iter
        
        % ============== ALMS iteration ===============
        close all; DispFilterSize=0; DispFilterStd=0; StrainFilterSize=0; StrainFilterStd=0; tol=1e-4;
        [UIter,FIter,alphaIter,normOfWIter,...
            coordinatesFEMIter,elementsFEMIter,eleGenerationIter,irregularEdgeIter,dirichletIter,neumannIter,...
            thetaDorfler,TimeICGN,EstimateTime,MarkTime,RefineTime,...
            rhoLastIterVector,rhoLastIterVectortemp1,rhoLastIterVectortemp2,EnrHAndTipEleIndex,EnrTipEleIndex] ...
            = funGBMSIterSq(UIter,FIter,coordinatesFEMIter,elementsFEMIter,winstepsize,ClusterNo,...
            eleGenerationIter,irregularEdgeIter,dirichletIter,neumannIter,LevelNo,ImgNormalized{1},ImgNormalized{2},winsize,...
            CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot,EnrHAndTipEleIndex,EnrTipEleIndex,...
            DispFilterSize,DispFilterStd,StrainFilterSize,StrainFilterStd,tol,alphaIter);
        
        LevelNo = LevelNo + 1; 
        FIter = funSmoothStrain_adapt(FIter,coordinatesFEMIter,elementsFEMIter,winstepsize,LevelNo,StrainFilterSize,StrainFilterStd);
  
        % ============== Assign values to variables ===============
        solnIter = struct('UIter',UIter,'FIter',FIter,'alphaIter',alphaIter,'normOfWIter',normOfWIter, ...
            'coordinatesFEMIter',coordinatesFEMIter,'elementsFEMIter',elementsFEMIter, ...
            'eleGenerationIter',eleGenerationIter,'dirichletIter',dirichletIter,'neumannIter',neumannIter, ...
            'TimeICGNIter',TimeICGN,'EstimateTimeLastIter',EstimateTime,'MarkTimeLastIter',MarkTime,'RefineTimeLastIter',RefineTime,...
            'rhoLastIterVector',rhoLastIterVector,'rhoLastIterVectortemp1',rhoLastIterVectortemp1,'rhoLastIterVectortemp2',rhoLastIterVectortemp2,'thetaDorflerLastIter',thetaDorfler,...
            'EnrHAndTipEleIndex',EnrHAndTipEleIndex,'EnrTipEleIndex',EnrTipEleIndex);
        
        soln{LevelNo} = solnIter;
        
        % ====== Plot refined mesh ======
        figure; show([], elementsFEMIter,coordinatesFEMIter,UIter(1:2:end));
        title('x-Disp U','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
        figure; show([], elementsFEMIter,coordinatesFEMIter,UIter(2:2:end));
        title('y-Disp V','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
        figure; show([], elementsFEMIter,coordinatesFEMIter,FIter(1:4:end));
        title('F11','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
        figure; show([], elementsFEMIter,coordinatesFEMIter,FIter(4:4:end));
        title('F22','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
        figure; show([], elementsFEMIter,coordinatesFEMIter,0.5*(FIter(2:4:end)+FIter(3:4:end)));
        title('0.5*(F12+F21)','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
        figure; show([], elementsFEMIter,coordinatesFEMIter,0.5*(FIter(2:4:end)-FIter(3:4:end)));
        title('0.5*(F21-F12)','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
        
        disp(['Level ',num2str(LevelNo),' is done. Press any key to keep refining mesh.']);
          
    end
    
     
   ResultSolnSeq{ImgSeqNum-1}=soln; %clear soln;
     
    
end




%%

DispFilterSize=0; DispFilterStd=0; StrainFilterSize=0; StrainFilterStd=0; tol=1e-4;tol2=1e-4;
 
LevelNo=5; ClusterNo=8;
solnIter = ResultSolnSeq{1}{LevelNo};
UIter = solnIter.UIter; FIter = solnIter.FIter; alphaIter = solnIter.alphaIter;     
coordinatesFEMIter = solnIter.coordinatesFEMIter; elementsFEMIter = solnIter.elementsFEMIter;
% irregularEdgeIter = solnIter.irregularEdgeIter;
eleGenerationIter = solnIter.eleGenerationIter; dirichletIter = solnIter.dirichletIter; neumannIter = solnIter.neumannIter;
EnrHAndTipEleIndex = solnIter.EnrHAndTipEleIndex; EnrTipEleIndex = solnIter.EnrTipEleIndex;
rhoLastIterVectorItertemp1 = solnIter.rhoLastIterVectortemp1; sum(rhoLastIterVectorItertemp1)


