% Matrix FDM 2D

function Lh = Matrix_FDM_2D_MorphogenGradient(N,h,k,gam)
%% Matrix_FDM_2D_MorphogenGradient
% This code calculates the matrix of the Laplacian operator for the
% Morphogen Gradient PDE problem studied in MCBU 2015. Check Vargas paper
% Vargas at. al. paper and Phay&Brian Notes to review the PDE.
% 
% cmd up to comment multiple line; cmd t to remove
clc;
% Parameters:
% k = 1;
% -------------------------------------------------------------------------
% Mesh Grid
% N = 3; % Number of points on the interval
% xm = -0.1;
% lm = 1;
% z1 = 0;
% zN = 1;
% gam = 1;
% hx = (lm-xm)/(N-1);
% hz = (zN-z1)/(N-1);
% h = max(hx,hz);
% k = 1;
%  Creating the sparse matrix that, for N = 3, looks like:
%          A   I   0
%    M  =  I   A   I
%          0   I   A
% 
% with I the identity matrix 3x3, 0 the zero matrix 3x3 and A the matrix
%          -4   1   0 
%    A  =   1  -4   1
%           0   1  -4

%  This matrix is the one you use when your BC are Dirichlet Boundary
%  conditions. If it is not the case, you modify the boundary nodes
%  accordingly with your BC.
I = speye(N);
e = ones(N,1);
A = spdiags([-e,4*e,-e],[-1,0,1],N,N);
J = spdiags([-e,-e],[-1,1],N,N);
Lh = kron(I,A) + kron(J,I);

% Modifying the boundary nodes. Recall that in our problem we have the 
% following boundary conditionsqQ: $$ \frac{\partial u}{\partial z} = 0 $$ at 
% z equals to 0, $$ \frac{\partial u}{\partial x} = 0 $$ at x = -xm, $$ u =
% 0 $$ at x = lm and $$ \frac{\partial u}{\partial z} + k*\gamma * u = 0 $$
% at z = 1.

%  For the corner point (-xm,0), so i =1, j = 1;
i = 1;
j = 1;
m = i + N*(j-1);
Lh(m,m+1) = -2;
Lh(m,m+N) = -2;
% For the boundary nodes (-xm,z), for z in (0,1), so i = 1, j = 2:N-1

i = 1;
for j = 2:N-1
    m = i + N*(j-1);
    Lh(m,m+1) = -2;
end
% For the corner point (-xm,1), so i = 1, j = N
i = 1;
j = N;
m = i + N*(j-1);
Lh(m,m) = 4+2*h*k*gam; 
Lh(m,m+1) = -2;
Lh(m,m-N) = -2;
% For the boundary nodes (x,1), so i = 2:N-1 and j = N
j = N;
for i = 2:N-1
    m = i + N*(j-1);
    Lh(m,m) = 4+2*h*k*gam;
    Lh(m,m-N) = -2;
end
% The boundary nodes (lm,z), for z in [0,1], so i = N and j = 1:N have the
% Dirichlet boundary condition u(i,j) = 0, so we not need to modify those
% nodes since the original matrix is for Dirichlet boundary conditions.
