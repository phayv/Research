
 function [bforcx,bforcy,sforcx,sforcy,sum2, mdr,theta]=forcedensity(m,dx,x,y,bendsti,ksurface,area0)
 
    fx=dx;
    x=x-mean(x);
    y=y-mean(y);
    
    dx=[x(1,2:m),x(1,1)]-x(1,1:m);
    dy=[y(1,2:m),y(1,1)]-y(1,1:m);
    
    ddx=[x(1,2:m),x(1,1)]-[x(1,m), x(1,1:m-1)];
    ddy=[y(1,2:m),y(1,1)]-[y(1,m), y(1,1:m-1)];
    
    
    kx=[x(1,2:m),x(1:1)]+[dx(1,m),dx(1,1:m-1)];
    ky=[y(1,2:m),y(1:1)]+[dy(1,m),dy(1,1:m-1)];
    
    dix=[x(1,2:m),x(1,1)]-kx;
    diy=[y(1,2:m),y(1,1)]-ky;
    
    dir=sign(-dy.*dix+dx.*diy);
    
    
    r=sqrt(dx.^2+dy.^2);
    mdr(1)=max(r);
    mdr(2)=min(r);
    r1=[r(1,m),r(1:m-1)];
    rr=sqrt(ddx.^2+ddy.^2);
    
    dtheta=(pi-real(acos((r.^2+r1.^2-rr.^2)./(2*r.*r1)))).*dir;
    
    rr1=r(1,1:m);
    rr2=-[r(1,m) r(1,1:m-1)];
    rr3=r(1,1:m)+[r(1,2:m) r(1,1)];
    
    dkap=dtheta./(r(1,1:m)-rr2)*2;
    dkap1=[dkap(1,2:m),dkap(1,1)];
    dkap2=[dkap(1,m),dkap(1,1:m-1),];
    dkap3=[dkap(1,3:m),dkap(1,1:2),];
    
    cf1=(rr2.*dkap1-rr1.*dkap2+(rr1-rr2).*dkap)./(rr2.*rr1).*(rr1.^2-rr3.^2)-...
        (rr3.*dkap1-rr1.*dkap3+(rr1-rr3).*dkap)./(rr3.*rr1).*(rr1.^2-rr2.^2);
    cf2=((rr1-rr2).*(rr1.^2-rr3.^2)-(rr1-rr3).*(rr1.^2-rr2.^2))/2;

    phi1=atan2(dy(1),dx(1));
    phi2=atan2(dy(m),dx(m));

    t0=(phi1+phi2)/2;

    stheta=(dtheta(1,1:m)+[dtheta(1,2:m),dtheta(1,1)])/2;

    theta(1,1:m+1)=cumsum([0,stheta(1,1:m)])+t0;

    rbn(1,1:m)=bendsti*ones(1,m);

    dkapss=cf1./cf2;

    y(m+1)=y(1);
    x(m+1)=x(1);
    sum2=sum(abs((x(1,2:m+1).*( -y(1,1:m))+x(1,1:m).*(y(1,2:m+1) )))/2);

    %farea=-ksurface*(sum2-area0)/fx/1000;    
    forc(1,1:m)=rbn(1,1:m).*(dkapss(1,1:m)+dkap(1,1:m).^3/2) ;% +farea;
    fforc=zeros(1,m);

    [bforcx,bforcy]=polortocartision(m,theta,forc,fforc);

    sforc=ksurface*(r-fx);

    sforcx=sforc.*(dx./r)+[sforc(1,m),sforc(1,1:m-1)].*[-dx(1,m),-dx(1,1:m-1)]./[r(1,m),r(1,1:m-1)];
    sforcy=sforc.*(dy./r)+[sforc(1,m),sforc(1,1:m-1)].*[-dy(1,m),-dy(1,1:m-1)]./[r(1,m),r(1,1:m-1)];

end
