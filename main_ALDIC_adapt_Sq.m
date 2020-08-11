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
try mex -O ba_interp2.cpp; catch; end % dbstop if error % % Old version codes.

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
fprintf('--- Extensions:  Additional parameters when dealing with DIC image sequence --- \n')
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ImgSeqNum = 2:length(ImgNormalized)
    
    % ============== Solve: Level1 uniform mesh ALDIC ==============
    [U0,ULoc,FLoc,USubpb2,FSubpb2,udual,vdual,ConvItPerEle,ALSolveStep,coordinatesFEM,elementsFEM,dirichlet,neumann,mu,beta,...
        ALSub1Time,ALSub2Time] = funALDIC(ImgNormalized,file_name,ImgSeqNum,gridxyROIRange,...
                                          winsize,winstepsize,ClusterNo,Subpb2FDOrFEM);
    dirichlet=[]; elementsFEM = [elementsFEM, zeros(size(elementsFEM))]; alpha = beta; eleGeneration = ones(size(elementsFEM,1),1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LevelNo = 1; UIter = full(USubpb2); FIter = full(FSubpb2); 
    udualIter = full(udual); vdualIter = full(vdual); 
    if CrackOrNot > 0, USubpb2tempIter = USubpb2temp;
    else, USubpb2tempIter = USubpb2; end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    % ============== Assign values to variables ===============
    ConvItPerEleIter = ConvItPerEle; betamuIter.betavar = beta; betamuIter.mu = mu; betamuIter.alphavar = alpha; 
    elementsFEMIter = elementsFEM; coordinatesFEMIter = coordinatesFEM; ALSolveStepIter = ALSolveStep;
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
    solnIter = struct('UIter',UIter,'FIter',FIter,'udualIter',udualIter,'vdualIter',vdualIter,...
        'USubpb2tempIter',USubpb2tempIter,'ConvItPerEleIter',ConvItPerEleIter,'betamuIter',betamuIter,...
        'coordinatesFEMIter',coordinatesFEMIter,'elementsFEMIter',elementsFEMIter,'dirichletIter',dirichletIter,'neumannIter',neumannIter,...
        'eleGenerationIter',eleGenerationIter, 'irregularEdgeIter', irregularEdgeIter, ... 
        'ALSub1TimeIter',ALSub1Time,'ALSub2TimeIter',ALSub2Time,'EstimateTimeLastIter',[0],'MarkTimeLastIter',[0],'RefineTimeLastIter',[0],...
        'rhoLastIterVector',[0],'rhoLastIterVectortemp1',[0],'rhoLastIterVectortemp2',[0],'thetaDorflerLastIter',[0],...
        'EnrHAndTipEleIndex',EnrHAndTipEleIndex,'EnrTipEleIndex',EnrTipEleIndex);
    
    soln{LevelNo} = solnIter;
     
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while (LevelNo<10) && (winstepsize/(2^(LevelNo-1))>1) && (mod((winstepsize/(2^(LevelNo-1))),2)==0)
      
        %% Adaptive-Refine - Iter
          
        % ============== ALMS iteration ===============
        close all; DispFilterSize=0; DispFilterStd=0; StrainFilterSize=0; StrainFilterStd=0; tol=1e-6;tol2=1e-4;
        [UIter,FIter,udualIter,vdualIter,USubpb2tempIter,ConvItPerEleIter,betamuIter,...
            coordinatesFEMIter,elementsFEMIter,eleGenerationIter,irregularEdgeIter,dirichletIter,neumannIter,...
            thetaDorfler,ALSub1Time,ALSub2Time,EstimateTime,MarkTime,RefineTime,...
            rhoLastIterVector,rhoLastIterVectortemp1,rhoLastIterVectortemp2,EnrHAndTipEleIndex,EnrTipEleIndex] ...
            = funALMSIterSq(UIter,FIter,udualIter,vdualIter,USubpb2tempIter,...
            coordinatesFEMIter,elementsFEMIter,winstepsize,ClusterNo,...
            eleGenerationIter,irregularEdgeIter,dirichletIter,neumannIter,LevelNo,ImgNormalized{1},ImgNormalized{2},winsize,...
            CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot,EnrHAndTipEleIndex,EnrTipEleIndex,...
            DispFilterSize,DispFilterStd,StrainFilterSize,StrainFilterStd,tol,tol2,betamuIter);
        
        % ============== Assign values to variables ===============
        solnIter = struct('UIter',UIter,'FIter',FIter,'udualIter',udualIter,'vdualIter',vdualIter,...
            'USubpb2tempIter',USubpb2tempIter,'ConvItPerEleIter',ConvItPerEleIter,'betamuIter',betamuIter,...
            'coordinatesFEMIter',coordinatesFEMIter,'elementsFEMIter',elementsFEMIter,'dirichletIter',dirichletIter,'neumannIter',neumannIter,...
        'eleGenerationIter',eleGenerationIter, 'irregularEdgeIter', irregularEdgeIter, ...  
            'ALSub1TimeIter',ALSub1Time,'ALSub2TimeIter',ALSub2Time,'EstimateTimeLastIter',EstimateTime,'MarkTimeLastIter',MarkTime,'RefineTimeLastIter',RefineTime,...
            'rhoLastIterVector',rhoLastIterVector,'rhoLastIterVectortemp1',rhoLastIterVectortemp1,'rhoLastIterVectortemp2',rhoLastIterVectortemp2,'thetaDorflerLastIter',thetaDorfler,...
            'EnrHAndTipEleIndex',EnrHAndTipEleIndex,'EnrTipEleIndex',EnrTipEleIndex);
        
        LevelNo = LevelNo + 1; soln{LevelNo} = solnIter;
         
        % ====== Plot refined mesh ======
%         figure; show([], elementsFEMIter(:,1:4),coordinatesFEMIter(:,1:2),UIter(1:2:end));
%         title('x-Disp U','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
%         figure; show([], elementsFEMIter(:,1:4),coordinatesFEMIter(:,1:2),UIter(2:2:end));
%         title('y-Disp V','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
%         figure; show([], elementsFEMIter(:,1:4),coordinatesFEMIter(:,1:2),FIter(1:4:end));
%         title('F11','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
%         figure; show([], elementsFEMIter(:,1:4),coordinatesFEMIter(:,1:2),FIter(4:4:end));
%         title('F22','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
%         figure; show([], elementsFEMIter(:,1:4),coordinatesFEMIter(:,1:2),0.5*(FIter(2:4:end)+FIter(3:4:end)));
%         title('0.5*(F12+F21)','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
%         figure; show([], elementsFEMIter(:,1:4),coordinatesFEMIter(:,1:2),0.5*(FIter(2:4:end)-FIter(3:4:end)));
%         title('0.5*(F21-F12)','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
      
        disp(['Level ',num2str(LevelNo),' is done. Press any key to keep refining mesh.']);
        pause;
        
    end
    
   
   ResultSolnSeq{ImgSeqNum-1}=soln; clear soln;
     
    
end

disp('--- Finish! ---')
%%


% figure,
DispFilterSize=0; DispFilterStd=0; StrainFilterSize=0; StrainFilterStd=0; tol=1e-4;tol2=1e-4;
 
LevelNo=4; ClusterNo=1;
solnIter = ResultSolnSeq{1}{LevelNo};
UIter = solnIter.UIter; FIter = solnIter.FIter; udualIter = solnIter.udualIter; vdualIter = solnIter.vdualIter;
USubpb2tempIter = solnIter.USubpb2tempIter; ConvItPerEleIter = solnIter.ConvItPerEleIter; betamuIter = solnIter.betamuIter;
coordinatesFEMIter = solnIter.coordinatesFEMIter; elementsFEMIter = solnIter.elementsFEMIter;
eleGenerationIter = solnIter.eleGenerationIter; dirichletIter = solnIter.dirichletIter; neumannIter = solnIter.neumannIter;
EnrHAndTipEleIndex = solnIter.EnrHAndTipEleIndex; EnrTipEleIndex = solnIter.EnrTipEleIndex;
rhoLastIterVectorItertemp1 = solnIter.rhoLastIterVectortemp1; sum(rhoLastIterVectorItertemp1)



% solnIter = ResultSolnSeq{1}{3};
%         UIter = solnIter.UIter; FIter = solnIter.FIter; udualIter = solnIter.udualIter; vdualIter = solnIter.vdualIter;
%         USubpb2tempIter = solnIter.USubpb2tempIter; ConvItPerEleIter = solnIter.ConvItPerEleIter; betamuIter = solnIter.betamuIter;
%         coordinatesFEMIter = solnIter.coordinatesFEMIter; elementsFEMIter = solnIter.elementsFEMIter; dirichletIter = solnIter.dirichletIter; neumannIter = solnIter.neumannIter;
%         eleGenerationIter = solnIter.eleGenerationIter; irregularEdgeIter = solnIter.irregularEdgeIter; 
%         EnrHAndTipEleIndex = solnIter.EnrHAndTipEleIndex; EnrTipEleIndex = solnIter.EnrTipEleIndex;
%          
% %close all;
% %Plotdisp_show(UIter,coordinatesFEMIter,elementsFEMIter(:,1:4));
% %    Plotstrain_show(FIter,coordinatesFEMIter,elementsFEMIter(:,1:4));
%     
% [row,col] = find(coordinatesFEMIter(:,2)==263);
%  coordinatesFEMtemp = coordinatesFEMIter(row,:);
%  F11Itertemp = FIter(4*row-3);
%  %UyItertemp = UIter(2*row);
% hold on;  plot(coordinatesFEMtemp(:,1),F11Itertemp,'k.')

 
