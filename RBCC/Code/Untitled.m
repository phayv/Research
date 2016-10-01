% THIS CODE COMPUTES THE SOLUTION TO THE clusetered RBCs in the capillary.

clear all;
%close all;
initialdata1=[128  1.*10^(-6)  1.*10^(-1)  1.*10^(-4)];
%            1   2     3    4
initialdata2=[.40  5*10^(-2)   1.8*10^(-19)   1*10^4];
%            1   2     3    4
initialdata3=[6 1  sqrt(2*1.38*10^(-23)*300*1.4*10^(-3)) 0];
%             1   2  3    4    5
initialdata4=[10^3  0  1  0];
%               1   2  3  4 
initialdata5=[1.6*10^(-3)  10^(-5)/2  0  11];
%               1          2        3  4
initialdata6=[5  3];
%             1  2

if initialdata3(4)==1
%     rng(initialdata5(4))
end

ngrid=round(initialdata1(1));      %the number of the grid points
dt=initialdata1(2);                %dt
T=initialdata1(3);                 %T
outpt=initialdata1(4);             %out put data groups

shortax=initialdata2(1);           %shape parameter
velocity=initialdata2(2);          %velocity
bendsti=initialdata2(3);           %bending stiffness
ksurface=initialdata2(4);          %elastic stiffness

Nvesicle=initialdata3(1);           %length of the pipe
Highth=initialdata3(2);            %highth of the pipe
Noi=initialdata3(3);               % noise
NOISEON=initialdata3(4);           % turn on noise or not

rho=initialdata4(1);               % mass density
INERTIAON=initialdata4(2);         %turn on inertia effect or not
Lam=initialdata4(3);               %viscosity ratio                
ONOFF=initialdata4(4);             %turn on or off viscosity ratio
                                   %  if ONOFF=1, viscosity is on, if==0
                                   %  off
mu=initialdata5(1);                %viscosity   
r0=initialdata5(2);                %typical radius of the vesicle
ktheta=initialdata5(3);		       %the initial tilted angle of the membrane

pressure=initialdata6(1);          % pressure 


shear=velocity;              
ktime=0;
m=ngrid;
Nv=Nvesicle;                


tcomp=dt*ktime;
nstep=round(T/dt);                %time steps
outpt=round(outpt/dt);            %steps for output
ksurface=ksurface/r0;

% %set up the intitial shape
r=1;
%[x,y,mr,sl]=equalinitial(m,r);
[x,y,mr,sl]=equalin(m,r);
 areao=areasum(m,x,y);
%x(m+1)=x(1);
%y(m+1)=y(1);
 
 R0=sqrt(areao/pi);
 delta=sl/R0-2*pi;


[x,y]=putintopipe(x,y,Nv);
%load('in6.mat')
Nv=6;
 [x,y,Nx,Ny,mr,dx]=scaleofpipe(x,y,mr,r0);
 areao=areasum(m,x,y);
Nx=round(Nx*2/3)-9;
 

sl=m*mr;
 
        xx(:,1)=x';
        yy(:,1)=y';
        are(:,1)=areao*ones(Nv,1);
         for i=1:Nv
            msr((i-1)*2+1:i*2,1)=[mr,mr]';               
         end             

WP=formWp(Nx,Ny,dx,mu);
