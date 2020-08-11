% Solving Ax = b by conjugate gradients.
% Reference:
% [1] https://www.mathworks.com/help/pdf_doc/otherdocs/simax.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = cgsolve (A,b,tol)
% Solve A*x = b by the conjugate gradient method.
% Iterate until norm(A*x-b) / norm(b) <= tol.
x = zeros(size(b));
r = b;
rtr = r'*r;
p = zeros(size(b));
beta = 0;

IterNo = 1; IterNoMax = 500;
while ( norm(r) > tol * norm(b) ) && ( IterNo < IterNoMax )
   
    p = r + beta * p;
    Ap = A * p;
    alpha = rtr / ( p' * Ap );
    x = x + alpha * p;
    r = r - alpha * Ap;
    rtrold = rtr;
    rtr = r'*r;
    beta = rtr / rtrold;
    
    IterNo = IterNo + 1;
    
end
