% ==============================================
% function funICGN 
% ==============================================
function [U,F,HtempPar,LocalTime,ConvItPerEle] = LocalICGN(U0,coordinatesFEM,...
    Df,imgfNormalizedbc,imggNormalizedbc,winsizeList,winstepsize,tol,ICGNmethod,ClusterNo)

ALSolveStep = 1;
temp = zeros(size(coordinatesFEM,1),1); UtempPar = temp; VtempPar = temp; % UtempPar = UPar{1}; VtempPar = UPar{2};
F11tempPar = temp; F21tempPar = temp; F12tempPar = temp; F22tempPar = temp;% F11tempPar = FPar{1}; F21tempPar = FPar{2}; F12tempPar = FPar{3}; F22tempPar = FPar{4};
HtempPar = zeros(size(coordinatesFEM,1),21);
ConvItPerEle = zeros(size(coordinatesFEM,1),1);
if length(winsizeList(:)) == 1, winsizeList = winsizeList*ones(size(coordinatesFEM,1),2);end

disp(['--- Start step ',num2str(ALSolveStep),' subproblem 1 ---']); tic;
% -------- How to change parallel pools ---------
% myCluster = parcluster('local');
% myCluster.NumWorkers = 4;  % 'Modified' property now TRUELevelNo
% saveProfile(myCluster);    % 'local' profile now updated,
%                            % 'Modified' property now FALSE
% -------------- Or we can do this --------------
% Go to the Parallel menu, then select Manage Cluster Profiles.
% Select the "local" profile, and change NumWorkers to 4.
% -----------------------------------------------
switch ClusterNo
case 0 || 1
    h = waitbar(0,'Please wait for Subproblem 1 IC-GN iterations!'); tic;
     
    for tempj = 1:size(coordinatesFEM,1)  % tempj is the element index
         try 
            winsize = winsizeList(tempj,:);
            x0temp = coordinatesFEM(tempj,1); y0temp = coordinatesFEM(tempj,2);  
            [Utemp, Ftemp, ConvItPerEle(tempj), HtempPar(tempj,:)] = funICGN(U0(2*tempj-1:2*tempj), ...
                    x0temp,y0temp,Df,imgfNormalizedbc,imggNormalizedbc,winsize,tol,ICGNmethod);
            %disp(['ele ',num2str(tempj),' converge step is ',num2str(ConvItPerEle(tempj)),' (>0-converged; 0-unconverged)']);
            % ------ Store solved deformation gradients ------
            UtempPar(tempj) = Utemp(1); VtempPar(tempj) = Utemp(2); 
            F11tempPar(tempj) = Ftemp(1); F21tempPar(tempj) = Ftemp(2); F12tempPar(tempj) = Ftemp(3); F22tempPar(tempj) = Ftemp(4);
            waitbar(tempj/(size(coordinatesFEM,1)));
         catch
             UtempPar(tempj) = nan; VtempPar(tempj) = nan;
             F11tempPar(tempj) = nan; F21tempPar(tempj) = nan; F12tempPar(tempj) = nan; F22tempPar(tempj) = nan;
           waitbar(tempj/(size(coordinatesFEM,1)));
        end
    end
    close(h); LocalTime = toc;
    
otherwise
    % Start parallel computing
    % ****** This step needs to be careful: may be out of memory ******
    tic;
    hbar = parfor_progressbar(size(coordinatesFEM,1),'Please wait for Subproblem 1 IC-GN iterations!');
    parfor tempj = 1:size(coordinatesFEM,1)  % tempj is the element index
        try 
            winsize = winsizeList(tempj,:);
            x0temp = coordinatesFEM(tempj,1); y0temp = coordinatesFEM(tempj,2);  
            [Utemp, Ftemp, ConvItPerEle(tempj), HtempPar(tempj,:)] = funICGN(U0(2*tempj-1:2*tempj), ...
                    x0temp,y0temp,Df,imgfNormalizedbc,imggNormalizedbc,winsize,tol,ICGNmethod);
            %disp(['ele ',num2str(tempj),' converge step is ',num2str(ConvItPerEle(tempj)),' (>0-converged; 0-unconverged)']);
            % ------ Store solved deformation gradients ------
            UtempPar(tempj) = Utemp(1); VtempPar(tempj) = Utemp(2); 
            F11tempPar(tempj) = Ftemp(1); F21tempPar(tempj) = Ftemp(2); F12tempPar(tempj) = Ftemp(3); F22tempPar(tempj) = Ftemp(4);
            hbar.iterate(1);
        catch
            UtempPar(tempj) = nan; VtempPar(tempj) = nan;
            F11tempPar(tempj) = nan; F21tempPar(tempj) = nan; F12tempPar(tempj) = nan; F22tempPar(tempj) = nan;
            hbar.iterate(1);
        end
    end
     
    close(hbar); 
    LocalTime = toc;
    
end

U = U0; U(1:2:end) = UtempPar; U(2:2:end) = VtempPar;
F = zeros(4*size(coordinatesFEM,1),1); F(1:4:end) = F11tempPar; F(2:4:end) = F21tempPar; F(3:4:end) = F12tempPar; F(4:4:end) = F22tempPar; 

% ------ Clear bad points for Local DIC ------
% find bad points after Local Subset ICGN
[row1,~] = find(ConvItPerEle(:)==0); 
[row2,~] = find(ConvItPerEle(:)>99);
row = unique(union(row1,row2));
disp(['Local step bad subsets total # is: ', num2str(length(row))]);
U(2*row-1) = NaN; U(2*row) = NaN; F(4*row-3) = NaN; F(4*row-2) = NaN;F(4*row-1) = NaN;F(4*row) = NaN;
% Plotdisp_show(full(USubpb1),elementsFEM(:,1:4) ,coordinatesFEM );
% Plotuv(full(USubpb1),x0,y0);
% ------ inpaint nans using gridfit ------
Coordxnodes = unique(coordinatesFEM(:,1)); Coordynodes = unique(coordinatesFEM(:,2));
nanindex = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);
[u1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); u1temp = u1temp';
[v1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex),Coordxnodes,Coordynodes,'regularizer','springs'); v1temp = v1temp';
[F11temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3),Coordxnodes,Coordynodes,'regularizer','springs'); F11temp = F11temp';
[F21temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2),Coordxnodes,Coordynodes,'regularizer','springs'); F21temp = F21temp';
[F12temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); F12temp = F12temp';
[F22temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0),Coordxnodes,Coordynodes,'regularizer','springs'); F22temp = F22temp';
  
% Add remove outliers median-test based on: 
% u1=cell(2,1); u1{1}=u1temp; u1{2}=v1temp;
% [u2] = removeOutliersMedian(u1,4); u2temp=u2{1}; v2temp=u2{2};
for tempi = 1:size(coordinatesFEM,1)
    [row1,col1] = find(Coordxnodes==coordinatesFEM(tempi,1));
    [row2,col2] = find(Coordynodes==coordinatesFEM(tempi,2));
    U(2*tempi-1) = u1temp(row1,row2);
    U(2*tempi)   = v1temp(row1,row2);
    F(4*tempi-3) = F11temp(row1,row2); F(4*tempi-2) = F21temp(row1,row2);
    F(4*tempi-1) = F12temp(row1,row2); F(4*tempi) = F22temp(row1,row2);
end

disp(['--- Done! step ',num2str(ALSolveStep),' subproblem 1 ---']);





