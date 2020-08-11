% ==================================
% function compute strain
% ----------------------------------
% switch MethodToComputeStrain
%   case 0: direct solved results;
%   case 1: central finite difference;
%   case 2: plane fitting method.
% ==================================

switch MethodToComputeStrain
    case 0
        FStrain = FLocal;  % JY!!! change intrisic coords to world coords
        FStraintemp = FStrain;
        Rad = 0;
        
        temp = 1:1:size(coordinatesFEM,1); temp = temp';
        temp = reshape(temp,M,N); temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad);
        temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad),1);
        
    case 1
        % Compute strain method I: Use Finite difference operator or FEM solver
        if Subpb2FDOrFEM == 1 %FD
			D = funDerivativeOp(M,N,winstepsize); % D = sparse(4*M*N, 2*M*N);
			FStrain = D*reshape(ULocal,length(ULocal),1); % JY!!! change intrisic coords to world coords
        else %FEM
			GaussPtOrder = 2; [FSubpb2,~,~] = funGlobal_NodalStrainAvg(coordinatesFEM,elementsFEM,USubpb2,GaussPtOrder);
        end
        % Plotstrain(FLocal,x,y,f,g);
        % End of computing strain method I.
        
        Rad = 1; % FilterSizeInput-2;
        % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
        temp = 1:1:size(coordinatesFEM,1); temp = temp';
        temp = reshape(temp,M,N); temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad);
        temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad),1);
        
        temp3 = zeros(4*(M-2*Rad)*(N-2*Rad),1);
        for i = 1:(M-2*Rad)*(N-2*Rad)
            temp3(4*i-3:4*i) = [4*temp2(i)-3; 4*temp2(i)-2; 4*temp2(i)-1; 4*temp2(i)];
        end
        FStraintemp = FStrain(temp3);
        
    case 2
        D = funDerivativeOp(M,N,winstepsize); % D = sparse(4*M*N, 2*M*N);
        FStrain = D*reshape(ULocal,length(ULocal),1); % JY!!! change intrisic coords to world coords
        
        % Compute strain method II: Use Plane Fitting method
        prompt = 'What is your half window size: ';
        Rad = input(prompt);
        
        [Uytemp, Uxtemp, UNewtemp] = PlaneFit2(reshape(ULocal(1:2:end),M,N), winstepsize, winstepsize,Rad);
        % figure; mesh(x((Rad+1):M-Rad,(Rad+1):N-Rad),y((Rad+1):M-Rad,(Rad+1):N-Rad),Uytemp((Rad+1):M-Rad,(Rad+1):N-Rad))
        % figure; mesh(x((Rad+1):M-Rad,(Rad+1):N-Rad),y((Rad+1):M-Rad,(Rad+1):N-Rad),Uxtemp((Rad+1):M-Rad,(Rad+1):N-Rad))
        % figure; mesh(x((Rad+1):M-Rad,(Rad+1):N-Rad),y((Rad+1):M-Rad,(Rad+1):N-Rad),UNewtemp((Rad+1):M-Rad,(Rad+1):N-Rad))
        [Vytemp, Vxtemp, VNewtemp] = PlaneFit2(reshape(ULocal(2:2:end),M,N), winstepsize, winstepsize,Rad);
        % figure; mesh(x((Rad+1):M-Rad,(Rad+1):N-Rad),y((Rad+1):M-Rad,(Rad+1):N-Rad),Vytemp((Rad+1):M-Rad,(Rad+1):N-Rad))
        % figure; mesh(x((Rad+1):M-Rad,(Rad+1):N-Rad),y((Rad+1):M-Rad,(Rad+1):N-Rad),Vxtemp((Rad+1):M-Rad,(Rad+1):N-Rad))
        % figure; mesh(x((Rad+1):M-Rad,(Rad+1):N-Rad),y((Rad+1):M-Rad,(Rad+1):N-Rad),VNewtemp((Rad+1):M-Rad,(Rad+1):N-Rad))
        
        FStraintemp = zeros(4*(M-2*Rad)*(N-2*Rad),1);
        FStraintemp(1:4:end) = reshape(Uxtemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        FStraintemp(2:4:end) = reshape(Vxtemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        FStraintemp(3:4:end) = reshape(Uytemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        FStraintemp(4:4:end) = reshape(Vytemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        
        % FStraintemp = -FStraintemp; % JY!!! change intrisic coords to world coords
        
        % ULocaltemp = zeros(2*(M-2*Rad)*(N-2*Rad),1);
        % ULocaltemp(1:2:end) = reshape(UNewtemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        % ULocaltemp(2:2:end) = reshape(VNewtemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        
        % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
        temp = 1:1:size(coordinatesFEM,1); temp = temp';
        temp = reshape(temp,M,N); temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad);
        temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad),1);
        
        % Plotstrain0(FStraintemp,x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad),f,g);
        
        temp3 = zeros(4*(M-2*Rad)*(N-2*Rad),1);
        for i = 1:(M-2*Rad)*(N-2*Rad)
            temp3(4*i-3:4*i) = [4*temp2(i)-3; 4*temp2(i)-2; 4*temp2(i)-1; 4*temp2(i)];
        end
        
        FStrain(temp3) = FStraintemp;
        
    otherwise
        disp('Wrong Input to compute strain field!')
        
end
