function [rhs,areas,mdrs,thetas]=rhsformPipe(m,dx,mr,x,y,Nx,Ny,Nv,noos,bendsti,ksurface,areao,pressure)
    press=pressure;
    Mn=length(x)/Nv;
    rhs=zeros(1,3*Nx*Ny);
 
    
    for i=1:Nv
        if i>1
            press=0;
        end
        
        noosI=noos(:,(i-1)*Mn+1:i*Mn);
        xI=x(1,(i-1)*Mn+1:i*Mn);
        yI=y(1,(i-1)*Mn+1:i*Mn);
        [bforcx,bforcy,sforcx,sforcy, sum2,mdr,theta]=forcedensity(m,mr,xI,yI,bendsti,ksurface,areao);
        
        rhs=rhs+rhsformp(m,Nx,Ny,dx,noosI,(-bforcx-sforcx)/dx,(-bforcy-sforcy)/dx,press);
        areas(i)=sum2;
        mdrs(i,:)=mdr;
        thetas(i,:)=theta(1,:);

    end



end

