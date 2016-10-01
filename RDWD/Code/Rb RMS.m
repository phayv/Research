% Rb RMS %
function Rb_RootsMeanSquare()

%--Boundary and Gridpoint Parameter------------------------%

N = 60; % Number of points on the interval
xm = -0.1;
lm = 1;
z1 = 0;
zN = 1;
hx = (lm-xm)/(N-1);
hz = (zN-z1)/(N-1);

%--Rb parameters-------------------------------------------%

Rb = 0; %initial guess
TolRb = 1e-5;
tol = 1e-5;
Rb_old = 1;
MaxIt = 100;
lm = 1;
zi = 1;
alpha = 1;
RB = zeros(MaxIt,1);
RB(1,1) = Rb;
[x,z] = meshgrid(xm:hx:lm, zN:-hz:z1);
x_old = x;
[~,n1] = size(x);
x = x >= 0;
[m,n2] = size(x);
n = n1 - n2;
x = reshape(x,m*n2,1);
%--------------------------------------------------------------------%
k = 1;
a_wild = main(tol,k,1,N);
b_wild = a_wild./(alpha+zi*a_wild);
bh = b_wild(1);
% bl = 0;
counter = 0;
while abs(Rb - Rb_old) > TolRb && counter < MaxIt
    e = 2;   
    a_ectopic = main(tol,k,e,N);
    b_ectopic = a_ectopic./(alpha+zi*a_ectopic);
% b_ectopic = a_ectopic./(...+zi*a_ectopic)
    Rb_old = Rb;
    Rb = (1/bh)*sqrt((1/lm)*trapz(x,(b_ectopic(m*n+1:end,1) - b_wild(n*m+1)).^2));
    k = update(Rb);    
    counter = counter + 1;
    RB(counter+1,1) = Rb;
end
b_ect_plot = reshape(b_ectopic,N,N);
b_wild_plot = reshape(b_wild,N,N);
figure(1)
surf(x_old,z,b_ect_plot)
xlabel('x')
ylabel('z')
zlabel('b_e_c_t_o_p_i_c')
figure(2)
surf(x_old,z,b_wild_plot)
figure(3)
contour(x_old,z,b_ect_plot)

figure(4)
contourf(x_old,z,b_ect_plot)

RB(1:counter+1,1)

%--Recalculates the kappa parameter---------------------------------%
    function k = update(Rb_old)
        c = 1;
        k = 1/(1+c*Rb_old);
    end
end
