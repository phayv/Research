



clear all;
%close all;
initialdata1=[128  1.*10^(-6)  1.*10^(-1)  1.*10^(-4)];
%            1   2     3    4
initialdata2=[.40   5*10^(-2)   1.8*10^(-19)   1*10^4];
%            1   2     3    4
initialdata3=[6 1  sqrt(2*1.38*10^(-23)*300*1.4*10^(-3)) 0];
%             1   2  3    4    5
initialdata4=[10^3  0  1  0];
%               1   2  3  4 
initialdata5=[1.6*10^(-3)  10^(-5)/2  0  11];
%               1          2        3  4
initialdata6=[2  3];
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

Nvesicle=initialdata3(1);            %length of the pipe
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
%Nx=Nx-mod(Nx,2);
%x=x-0.02*max(x);

% plot(x,y,'*')
%     lsy=[dx:dx:dx*Ny]-dx*Ny/2-dx/2;
%     lsx=[dx:dx:dx*Nx]-dx*Nx/2-dx/2;
%     hold on;plot(lsx,0,'ro');
%     plot(0,lsy,'g<')
% Nx=Nx*2;
% Ny=Ny*2;
% dx=dx/2;
% Nv=4;
%plot(x,y,'.');axis image
sl=m*mr;
%Noi=Noi/sqrt(dt*dx^2);
%%% get the initial data
        xx(:,1)=x';
        yy(:,1)=y';
        are(:,1)=areao*ones(Nv,1);
         for i=1:Nv
            msr((i-1)*2+1:i*2,1)=[mr,mr]';               
         end             

%form the right hand side by the latice      
W=formWp(Nx,Ny,1,1);
%W=formWp(Nx,Ny,dx,mu);
%W=formWinkle(Nx,Ny,dx,mu);
WP=W;          

coe1=dx/mu;
coe2=dx^2/mu;

[x,y]=vesiclereassign(m,x,y,-1);
x=x+2*dx;


for ktime=1:nstep
    
        if max(x)>=(Nx/2-2)*dx;
            x=perioddet(Nv,x,Nx,dx);
        end
    %ktime
    noos=directperid(m,Nx,Ny,x,y,dx);
    [rhs,areas,mdrs,thetas]=rhsformPipe(m,dx,mr,x,y,Nx,Ny,Nv,noos,bendsti,ksurface,areao,pressure);
    vel= WP\(coe2*rhs)';

    vell=vel';           
    unn=vell(1:Nx*Ny);
    utt=vell(Nx*Ny+1:2*Nx*Ny);
    P=vell(2*Nx*Ny+1:3*Nx*Ny);                  
    %put unn into ux uy
    
    [uxs,uys]=deltaprimep(Nx,Ny,dx,noos,unn,utt);

    if ktime==1
        ux=uxs;
        uy=uys;        
    end
    
    x=x+uxs*dt*3/2-ux*dt/2;          
    y=y+uys*dt*3/2-uy*dt/2;
 
     ux=uxs;
     uy=uys;

     if mod(ktime,outpt)==0
         tcomp=dt*ktime
         np=ktime/outpt+1;
         xx(:,np)=x';
         yy(:,np)=y';
         tt(np)=tcomp;
         are(:,np)=areas;   
         for i=1:Nv
            msr(i,np)=mean(mdrs(i,:)');               
         end
         xx=real(xx);
         yy=real(yy);
         x=real(x);
         y=real(y);
         noos=real(noos);
     end              
 % plot(x,y,'r'); hold on;
  if  mod(ktime,10000)==0
      save('w6p-128-k18-co40000-mu16-p2-1dt6.mat',...
	'xx','yy','tt','are','msr','vel','dt','T',...
	'noos','dx','ksurface','velocity','Noi','Nx','Ny');
  end
end






