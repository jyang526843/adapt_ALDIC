 function A = funDerivativeOp(M,N,h)

A = sparse(4*M*N, 2*M*N);
A(1,1) = -2; A(1,3)= 2; 
A(3,1) = -2; A(3,2*(M+1)-1) = 2;
A(2,2) = -2; A(2,4) = 2;
A(4,2) = -2; A(4,2*(M+1)) = 2;

% Inside
for tempx = 2:(M-1)
    for tempy = 2:(N-1)
        index = tempx+M*(tempy-1); % Find the point position
        indexUp = tempx+M*tempy;
        indexDown = tempx+M*(tempy-2);
        indexLeft = tempx-1+M*(tempy-1);
        indexRight = tempx+1+M*(tempy-1);
        
        indexF11row = 4*index-3; 
        indexF21row = 4*index-2;
        indexF12row = 4*index-1;
        indexF22row = 4*index;
        
        indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(indexRight)-1;
        indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(indexUp)-1;
        indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(indexRight);
        indexF22col1 = 2*(indexDown); indexF22col2 = 2*(indexUp);
        
        A(indexF11row,indexF11col1) = -1; A(indexF11row,indexF11col2) = 1;
        A(indexF21row,indexF21col1) = -1; A(indexF21row,indexF21col2) = 1;
        A(indexF12row,indexF12col1) = -1; A(indexF12row,indexF12col2) = 1;
        A(indexF22row,indexF22col1) = -1; A(indexF22row,indexF22col2) = 1; 
    end
end

% Bottom line
for tempx = 2:(M-1)
    for tempy = 1
        index = tempx+M*(tempy-1); % Find the point position
        indexUp = tempx+M*tempy;
        % indexDown = tempx+M*(tempy-2);
        indexLeft = tempx-1+M*(tempy-1);
        indexRight = tempx+1+M*(tempy-1);
        
        indexF11row = 4*index-3; 
        indexF21row = 4*index-2;
        indexF12row = 4*index-1;
        indexF22row = 4*index;
        
        indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(indexRight)-1;
        indexF12col1 = 2*(index)-1; indexF12col2 = 2*(indexUp)-1;
        indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(indexRight);
        indexF22col1 = 2*(index); indexF22col2 = 2*(indexUp);
        
        A(indexF11row,indexF11col1) = -1; A(indexF11row,indexF11col2) = 1;
        A(indexF21row,indexF21col1) = -1; A(indexF21row,indexF21col2) = 1;
        A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
        A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2; 
        
    end
end

% Top
for tempx = 2:(M-1)
    for tempy = N
        index = tempx+M*(tempy-1); % Find the point position
        % indexUp = tempx+M*tempy;
        indexDown = tempx+M*(tempy-2);
        indexLeft = tempx-1+M*(tempy-1);
        indexRight = tempx+1+M*(tempy-1);
        
        indexF11row = 4*index-3; 
        indexF21row = 4*index-2;
        indexF12row = 4*index-1;
        indexF22row = 4*index;
        
        indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(indexRight)-1;
        indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(index)-1;
        indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(indexRight);
        indexF22col1 = 2*(indexDown); indexF22col2 = 2*(index);
        
        A(indexF11row,indexF11col1) = -1; A(indexF11row,indexF11col2) = 1;
        A(indexF21row,indexF21col1) = -1; A(indexF21row,indexF21col2) = 1;
        A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
        A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2; 
    end
end

% Left
for tempx = 1
    for tempy = 2:(N-1)
        index = tempx+M*(tempy-1); % Find the point position
        indexUp = tempx+M*tempy;
        indexDown = tempx+M*(tempy-2);
        % indexLeft = tempx-1+M*(tempy-1);
        indexRight = tempx+1+M*(tempy-1);
        
        indexF11row = 4*index-3; 
        indexF21row = 4*index-2;
        indexF12row = 4*index-1;
        indexF22row = 4*index;
        
        indexF11col1 = 2*(index )-1; indexF11col2 = 2*(indexRight)-1;
        indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(indexUp)-1;
        indexF21col1 = 2*(index ); indexF21col2 = 2*(indexRight);
        indexF22col1 = 2*(indexDown); indexF22col2 = 2*(indexUp);
        
        A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
        A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
        A(indexF12row,indexF12col1) = -1; A(indexF12row,indexF12col2) = 1;
        A(indexF22row,indexF22col1) = -1; A(indexF22row,indexF22col2) = 1; 
    end
end

% Right
for tempx = M
    for tempy = 2:(N-1)
        index = tempx+M*(tempy-1); % Find the point position
        indexUp = tempx+M*tempy;
        indexDown = tempx+M*(tempy-2);
        indexLeft = tempx-1+M*(tempy-1);
        % indexRight = tempx+1+M*(tempy-1);
        
        indexF11row = 4*index-3; 
        indexF21row = 4*index-2;
        indexF12row = 4*index-1;
        indexF22row = 4*index;
        
        indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(index )-1;
        indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(indexUp)-1;
        indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(index );
        indexF22col1 = 2*(indexDown); indexF22col2 = 2*(indexUp);
        
        A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
        A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
        A(indexF12row,indexF12col1) = -1; A(indexF12row,indexF12col2) = 1;
        A(indexF22row,indexF22col1) = -1; A(indexF22row,indexF22col2) = 1; 
    end
end


% LeftBottom
for tempx = 1
    for tempy = 1
        index = tempx+M*(tempy-1); % Find the point position
        indexUp = tempx+M*tempy;
        % indexDown = tempx+M*(tempy-2);
        % indexLeft = tempx-1+M*(tempy-1);
        indexRight = tempx+1+M*(tempy-1);
        
        indexF11row = 4*index-3; 
        indexF21row = 4*index-2;
        indexF12row = 4*index-1;
        indexF22row = 4*index;
        
        indexF11col1 = 2*(index )-1; indexF11col2 = 2*(indexRight)-1;
        indexF12col1 = 2*(index )-1; indexF12col2 = 2*(indexUp)-1;
        indexF21col1 = 2*(index ); indexF21col2 = 2*(indexRight);
        indexF22col1 = 2*(index ); indexF22col2 = 2*(indexUp);
        
        A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
        A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
        A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
        A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2; 
    end
end

% RightBottom
for tempx = M
    for tempy = 1
        index = tempx+M*(tempy-1); % Find the point position
        indexUp = tempx+M*tempy;
        %indexDown = tempx+M*(tempy-2);
        indexLeft = tempx-1+M*(tempy-1);
        %indexRight = tempx+1+M*(tempy-1);
        
        indexF11row = 4*index-3; 
        indexF21row = 4*index-2;
        indexF12row = 4*index-1;
        indexF22row = 4*index;
        
        indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(index )-1;
        indexF12col1 = 2*(index )-1; indexF12col2 = 2*(indexUp)-1;
        indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(index );
        indexF22col1 = 2*(index ); indexF22col2 = 2*(indexUp);
        
        A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
        A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
        A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
        A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2; 
    end
end

% LeftTop
for tempx = 1
    for tempy = N
        index = tempx+M*(tempy-1); % Find the point position
        %indexUp = tempx+M*tempy;
        indexDown = tempx+M*(tempy-2);
        %indexLeft = tempx-1+M*(tempy-1);
        indexRight = tempx+1+M*(tempy-1);
        
        indexF11row = 4*index-3; 
        indexF21row = 4*index-2;
        indexF12row = 4*index-1;
        indexF22row = 4*index;
        
        indexF11col1 = 2*(index )-1; indexF11col2 = 2*(indexRight)-1;
        indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(index )-1;
        indexF21col1 = 2*(index ); indexF21col2 = 2*(indexRight);
        indexF22col1 = 2*(indexDown); indexF22col2 = 2*(index );
        
        A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
        A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
        A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
        A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2; 
    end
end

% RightTop
for tempx = M
    for tempy = N
        index = tempx+M*(tempy-1); % Find the point position
        %indexUp = tempx+M*tempy;
        indexDown = tempx+M*(tempy-2);
        indexLeft = tempx-1+M*(tempy-1);
        %indexRight = tempx+1+M*(tempy-1);
        
        indexF11row = 4*index-3; 
        indexF21row = 4*index-2;
        indexF12row = 4*index-1;
        indexF22row = 4*index;
        
        indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(index)-1;
        indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(index)-1;
        indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(index);
        indexF22col1 = 2*(indexDown); indexF22col2 = 2*(index);
        
        A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
        A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
        A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
        A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2; 
    end
end


A = A/(2*h);