function [x1,y1,mr,sl]=equalin(m,ynot)

for i=1:m
    theta(i)=2*pi*(i-1)/m;
    x(i)=ynot/2.8324*cos(theta(i))+ynot/3*2*cos(2*theta(i));
    y(i)=ynot*sin(theta(i))+ynot/3*0.6*sin(2*theta(i));
end

%plot(x,y,'.'); axis image
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


end