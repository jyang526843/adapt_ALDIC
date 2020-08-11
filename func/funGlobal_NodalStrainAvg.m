function [Strain] = funGlobal_NodalStrainAvg(coordinates,elements,U,GaussPtOrder)

FStrainAvgTimes = zeros(4*size(coordinates,1),1);
FStrain = zeros(4*size(coordinates,1),1);

for j = 1:size(elements,1) % j is the element index
    
    StrainWithinEachElementGausspoint = zeros(4,4);
    
    % ----- Find four corner points -------
    point1x = coordinates(elements(j,1),1);
    point1y = coordinates(elements(j,1),2);
    point2x = coordinates(elements(j,2),1);
    point2y = coordinates(elements(j,2),2);
    point3x = coordinates(elements(j,3),1);
    point3y = coordinates(elements(j,3),2);
    point4x = coordinates(elements(j,4),1);
    point4y = coordinates(elements(j,4),2);
    
    % ------ Find the element nodal indices ------
    temp = [2*elements(j,1)-1 2*elements(j,1) 2*elements(j,2)-1 2*elements(j,2)...
        2*elements(j,3)-1 2*elements(j,3) 2*elements(j,4)-1 2*elements(j,4)];
    
    % ------ Four Gauss points ------
    for tempi = 1:2
        for tempj = 1:2
            
            %ksi = ( 2*tempj-3)/sqrt(3); eta = (2*tempi-3)/sqrt(3);
            ksi = 2*tempj-3; eta = 2*tempi-3;
            if (tempi == 1) && (tempj == 1)
                ksi = -1/sqrt(3); eta = -1/sqrt(3);
            elseif (tempi == 1) && (tempj == 2)
                ksi = 1/sqrt(3); eta = -1/sqrt(3);
            elseif (tempi == 2) && (tempj == 1)
                ksi = 1/sqrt(3); eta = 1/sqrt(3);
            elseif (tempi == 2) && (tempj == 2)
                ksi = -1/sqrt(3); eta = 1/sqrt(3);
            end
            
            
            % ------ Calculate N ------
            N1 = (1-ksi)*(1-eta)*0.25;
            N2 = (1+ksi)*(1-eta)*0.25;
            N3 = (1+ksi)*(1+eta)*0.25;
            N4 = (1-ksi)*(1+eta)*0.25;
            N = [N1 0 N2 0 N3 0 N4 0;
                0 N1 0 N2 0 N3 0 N4];
            
            
            % ------ Build J matrix ------
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
            
            % ------ Compute DN matrix ------
            DN = [InvJ zeros(2,2); zeros(2,2) InvJ] * ...
                [funDN1x(ksi,eta) 0 funDN2x(ksi,eta) 0 funDN3x(ksi,eta) 0 funDN4x(ksi,eta) 0;
                funDN1y(ksi,eta) 0 funDN2y(ksi,eta) 0 funDN3y(ksi,eta) 0 funDN4y(ksi,eta) 0;
                0 funDN1x(ksi,eta) 0 funDN2x(ksi,eta) 0 funDN3x(ksi,eta) 0 funDN4x(ksi,eta);
                0 funDN1y(ksi,eta) 0 funDN2y(ksi,eta) 0 funDN3y(ksi,eta) 0 funDN4y(ksi,eta)];
            
            StrainWithinEachElementGausspoint(2*(tempi-1)+tempj,1:4) = DN*U(temp);
            
            
        end
    end
    
    MatrixExtrapolation = [1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3) -0.5;
        -0.5 1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3);
        1-0.5*sqrt(3) -0.5 1+0.5*sqrt(3) -0.5;
        -0.5 1-0.5*sqrt(3) -0.5 1+0.5*sqrt(3)];
    % ------ Nodal points strain extrapolation using Gauss points -----
    StrainExxWithinEachElementNodalpoint =  MatrixExtrapolation * [StrainWithinEachElementGausspoint(1,1);
        StrainWithinEachElementGausspoint(2,1);
        StrainWithinEachElementGausspoint(3,1);
        StrainWithinEachElementGausspoint(4,1)];
    
    StrainExyWithinEachElementNodalpoint = MatrixExtrapolation * [StrainWithinEachElementGausspoint(1,2);
        StrainWithinEachElementGausspoint(2,2);
        StrainWithinEachElementGausspoint(3,2);
        StrainWithinEachElementGausspoint(4,2)];
    
    StrainEyxWithinEachElementNodalpoint = MatrixExtrapolation * [StrainWithinEachElementGausspoint(1,3);
        StrainWithinEachElementGausspoint(2,3);
        StrainWithinEachElementGausspoint(3,3);
        StrainWithinEachElementGausspoint(4,3)];
    
    StrainEyyWithinEachElementNodalpoint = MatrixExtrapolation * [StrainWithinEachElementGausspoint(1,4);
        StrainWithinEachElementGausspoint(2,4);
        StrainWithinEachElementGausspoint(3,4);
        StrainWithinEachElementGausspoint(4,4)];
    
    StrainWithinEachElementGausspoint(1,1) = StrainExxWithinEachElementNodalpoint(1);
    StrainWithinEachElementGausspoint(2,1) = StrainExxWithinEachElementNodalpoint(2);
    StrainWithinEachElementGausspoint(3,1) = StrainExxWithinEachElementNodalpoint(3);
    StrainWithinEachElementGausspoint(4,1) = StrainExxWithinEachElementNodalpoint(4);
    
    StrainWithinEachElementGausspoint(1,2) = StrainExyWithinEachElementNodalpoint(1);
    StrainWithinEachElementGausspoint(2,2) = StrainExyWithinEachElementNodalpoint(2);
    StrainWithinEachElementGausspoint(3,2) = StrainExyWithinEachElementNodalpoint(3);
    StrainWithinEachElementGausspoint(4,2) = StrainExyWithinEachElementNodalpoint(4);
    
    StrainWithinEachElementGausspoint(1,3) = StrainEyxWithinEachElementNodalpoint(1);
    StrainWithinEachElementGausspoint(2,3) = StrainEyxWithinEachElementNodalpoint(2);
    StrainWithinEachElementGausspoint(3,3) = StrainEyxWithinEachElementNodalpoint(3);
    StrainWithinEachElementGausspoint(4,3) = StrainEyxWithinEachElementNodalpoint(4);
    
    StrainWithinEachElementGausspoint(1,4) = StrainEyyWithinEachElementNodalpoint(1);
    StrainWithinEachElementGausspoint(2,4) = StrainEyyWithinEachElementNodalpoint(2);
    StrainWithinEachElementGausspoint(3,4) = StrainEyyWithinEachElementNodalpoint(3);
    StrainWithinEachElementGausspoint(4,4) = StrainEyyWithinEachElementNodalpoint(4);
    
    
    % ------ Find the element nodal indices for strain ------
    tempStrainIndex = [4*elements(j,1)-3 4*elements(j,1)-2  4*elements(j,1)-1 4*elements(j,1)  ...
        4*elements(j,2)-3 4*elements(j,2)-2  4*elements(j,2)-1 4*elements(j,2) ...
        4*elements(j,3)-3 4*elements(j,3)-2  4*elements(j,3)-1 4*elements(j,3) ...
        4*elements(j,4)-3 4*elements(j,4)-2  4*elements(j,4)-1 4*elements(j,4)];
    
    
    FStrain(tempStrainIndex) = FStrain(tempStrainIndex)+[StrainWithinEachElementGausspoint(1,1);
        StrainWithinEachElementGausspoint(1,2); StrainWithinEachElementGausspoint(1,3);
        StrainWithinEachElementGausspoint(1,4); StrainWithinEachElementGausspoint(2,1);
        StrainWithinEachElementGausspoint(2,2); StrainWithinEachElementGausspoint(2,3);
        StrainWithinEachElementGausspoint(2,4); StrainWithinEachElementGausspoint(3,1);
        StrainWithinEachElementGausspoint(3,2); StrainWithinEachElementGausspoint(3,3);
        StrainWithinEachElementGausspoint(3,4); StrainWithinEachElementGausspoint(4,1);
        StrainWithinEachElementGausspoint(4,2); StrainWithinEachElementGausspoint(4,3);
        StrainWithinEachElementGausspoint(4,4)];
    
    FStrainAvgTimes(tempStrainIndex) = FStrainAvgTimes(tempStrainIndex) + ones(16,1);
    
    
    
    
end


Strain = FStrain./FStrainAvgTimes;
                       
            
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
