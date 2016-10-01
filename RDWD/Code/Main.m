% Main %
function u = main(tol,k,e,N)
%% Newton method

%----------------Modify Parameters-----------------------------------%

%Boundary and gridpoint parameters
% N = 30; % Number of points on the interval
xm = -0.1;
lm = 1;
z1 = 0;
zN = 1;
gam = 1;
hx = (lm-xm)/(N-1);
hz = (zN-z1)/(N-1);
h = max(hx,hz);

%Constant parameters, meshgrid point and max iterations
MaxIt = 30;
% k = 1;
gz = 1;
alpha = 1;
zi = 1;
vz = 1;
[x,z] = meshgrid(xm:hx:lm, zN:-hz:z1);
theta = 0.6;

%--Functions---------------------------------------------------------%

u = zeros(N^2,1); % vector of initial guess

Lh = Matrix_FDM_2D_MorphogenGradient(N,h,k,gam); %Lapacian Operator Matrix

f = @(u) - k^2*gz*u/(alpha + zi*u);
fp = @(u) - k^2*gz/(alpha + zi*u) + zi*k^2*gz*u/(alpha + zi*u)^2;
v = @(x) e*k^2*vz*heaviside(-x);
F = RHS_FDM_2D_MorphogenGradient(N,h,f,u,v,x);
Fp = derivativeRHS_FDM_2D_Morphogen(N,h,fp,u);
%--------------------------------------------------------------------%
% Newton's Method
 r0 = F;
 r = r0;
 counter = 0;
 while norm(r)/norm(r0) >= tol && counter < MaxIt
       r = F - Lh*u;
       Fp = spdiags(Fp,0,N^2,N^2);
       J = Lh - Fp;
       e = inexact(J,r,theta);
       u_old = u;
       u = update(u_old,e);
       counter = counter + 1;
 end
 
 u_plot = reshape(u,N,N);
 surf(x,z,u_plot)
 
% --Inexact Method---------------------------------------------------%
        function e = inexact(J,r,theta)
            e = zeros(N^2,1); 
            r_in = r;
            counter1 = 0;
            while norm(r_in) > theta*norm(r) && counter1 < MaxIt
                e = e + J\r;
                r_in = r -J*e;
                counter1 = counter1 + 1;
            end
        end

% --Updating a new u vector------------------------------------------%        
        function u = update(u_old,e)
            s = 1; % Simple bisection scheme
            r_old = norm(Lh*u_old - RHS_FDM_2D_MorphogenGradient(N,h,f,u_old,v,x));
            r = r_old + 1; % Error Tolerance
            counter2 = 0;
            while r > r_old && counter2 < MaxIt
                 u = u_old + s*e;
                 r = norm(Lh*u - F);
                 s = s/2;
                 counter2 = counter2 + 1;
            end
        end
end
