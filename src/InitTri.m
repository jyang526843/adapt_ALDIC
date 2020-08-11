%% Correct and Smooth initial displacements
function U0 =InitTri(u,v,Phi,uCen,vCen,PhiCen,M,N)


% ========= Delete bad intial points ==========
% In inpaint_nans part, I use Spring Model for u and v initial guesses
uInit = u; vInit = v; 
threshod = 0.3;
[row, col] = find(Phi<(1-threshod));
uInit(row,col) = NaN; vInit(row,col) = NaN;

uInit = inpaint_nans(uInit,4);
vInit = inpaint_nans(vInit,4);

try
[row, col] = find(PhiCen<(1-threshod));
uCen(row,col) = NaN; vCen(row,col) = NaN;

uCen = inpaint_nans(uCen,4);
vCen = inpaint_nans(vCen,4);
catch
end

% ========= Set initial values for U =========
U0 = zeros(2*(M*N+(M-1)*(N-1)),1);
uInit = uInit'; vInit = vInit'; % Here transpose bcz following give valus in column order
uCen = uCen';  vCen = vCen';
for tempi = 1:M*N  
    U0(2*tempi-1) = uInit(tempi);  
    U0(2*tempi)   = vInit(tempi);  
end

try
for tempi = 1:((M-1)*(N-1))
    U0(2*(M*N+tempi)-1) = uCen(tempi);
    U0(2*(M*N+tempi))   = vCen(tempi);
end
catch
end

disp('Finish setting up mesh and assigning initial value!')
  