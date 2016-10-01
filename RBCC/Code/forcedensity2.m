function [bforcx,bforcy,sforcx,sforcy,sum2, mdr,theta]=forcedensity2(m,dx,x,y,bendsti,ksurface,area0)
        %[bforcx,bforcy,sforcx,sforcy,sum2,mdr]=bendingforcesp(m,sl,x,y,bendsti,ksurface,area0)
    fx=dx;
    coef1=[0 : m/2  -m/2+1:-1]*sqrt(-1)* (2*pi)/m;
    coef12=coef1.^2;
    coef2=[1 1:m/2  -m/2+1:-1]*sqrt(-1)* (2*pi)/m;
    
    
    
    xk=fft(x,m);
    xk1=xk.*coef1 ;
    xk2=-xk.*coef12;
    
    xd1=real(ifft(xk1,m));
    xd2=real(ifft(xk2,m));
    
    yk=fft(y,m);
    yk1=yk.*coef1;
    yk2=-yk.*coef12;
    
    yd1=real(ifft(yk1,m));
    yd2=real(ifft(yk2,m));    
    
    sk=sqrt(xd1.^2+yd1.^2);
    
    dkap=(xd2.*yd1-xd1.*yd2)./sk.^3;
    
 
    dkapk=fft(dkap,m);
    dkapk=dkapk.*coef1;
    
    dkapk1=real(ifft(dkapk,m)); 
    
    dkaps1=dkapk1./sk;
    
    
    dkapk=fft(dkaps1,m);
    dkapk=dkapk.*coef1;
    dkapk2=real(ifft(dkapk,m)); 
    
    dkapss=dkapk2./sk;    
 
    thetak=dkap.*sk;
    
    thetak=fft(thetak,m)./coef2;
    
    theta=real(ifft(thetak,m));
    
    dte=2*pi/m;
    
    theta=theta+[0:m-1]*dte;
    
    dx=[x(1,2:m),x(1,1)]-x(1,1:m);
    dy=[y(1,2:m),y(1,1)]-y(1,1:m);
    r=sqrt(dx.^2+dy.^2);
    phi=atan2(dy,dx);
    
    for i=1:m-1
        if phi(i+1)-phi(i)<-pi
            phi(i+1:m)=phi(i+1:m)+2*pi;
        end
    end
    
    psi=(phi(1:m)+[phi(m)-2*pi phi(1:m-1)])/2;
    
    theta=theta+mean(psi-theta);
    
    
    rbn(1,1:m)=bendsti*ones(1,m);
    forc(1,1:m)=rbn(1,1:m).*(dkapss(1,1:m)+dkap(1,1:m).^3/2);
    fforc=zeros(1,m);    
    [bforcx,bforcy]=polortocartision(m,theta,forc,fforc);
    
    
    sd=(sk+[sk(2:m) sk(1)])/2;
    
    sforc=8*ksurface*(sd-fx);
    
    sforcx=sforc.*(dx./r)+[sforc(1,m),sforc(1,1:m-1)].*[-dx(1,m),-dx(1,1:m-1)]./[r(1,m),r(1,1:m-1)];
    sforcy=sforc.*(dy./r)+[sforc(1,m),sforc(1,1:m-1)].*[-dy(1,m),-dy(1,1:m-1)]./[r(1,m),r(1,1:m-1)];
    
 
    
    sum2= sum((sin(theta).*x(1,1:m)-cos(theta).*y(1,1:m)).*sk);
    
    mdr=sk;
    
    

end
