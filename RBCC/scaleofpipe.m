function [x0,y0,Nx,Ny,mr1,dx]=scaleofpipe(x,y,mr,r0)

    XL=max(abs(x));
    Yl=max(abs(y));
    r=r0/2;
    coef=r/Yl;
    
    x0=x*coef;
    y0=y*coef;
    mr1=mr*coef;
    
    Ny=ceil(2.4*r/mr1);
    if mod(Ny,2)==1
        Ny=Ny+1;
    end
    dx=2.4*r/(Ny-1);
    
    Nx=ceil(3*XL*coef/dx);
    
    if mod(Nx,2)==1
        Nx=Nx+1;
    end
    

end