function [x1,y1,mr,sl]=equalinitial(m,r)

for i=1:m
    theta(i)=2*pi*(i-1)/m;
    x(i)=(1+0.8*cos(pi+2*theta(i)))*cos(theta(i));
    y(i)=(2+0.4*cos(pi+2*theta(i)))*sin(theta(i));
end

Ly=max(y);
x=r*x/Ly;
y=r*y/Ly;
x1=x;
y1=y;
ppx = csape([0:m],[x,x(1)],'periodic');
ppy = csape([0:m],[y,y(1)],'periodic');

nm=[0:m-1];
for k=1:200    
    x1= ppval(ppx,nm);
    y1= ppval(ppy,nm);    
    %hold on; plot(x1,y1,'ro'); axis image
    dx=[x1(1,2:m),x1(1,1)]-x1(1,1:m);
    dy=[y1(1,2:m),y1(1,1)]-y1(1,1:m);
    dr=sqrt(dx.^2+dy.^2);
    mr=mean(dr);
    essor(k)=max(abs(dr-mr)/mr);
    if essor(k)<10^(-12)
        k
        break
    end    
    kr=cumsum([0,dr]);    
    %plot(x1,y1,'ro'); hold on
    pp1 = csape(kr,[nm,nm(1)+m],'periodic');
    nm = ppval(pp1,[0:m-1]*mr);
    
end

sl=m*mr;
%plot(x1,y1,'*');axis image
end