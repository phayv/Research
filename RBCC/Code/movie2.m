function t=movie2(x,y,dx,Nx,Ny,t,kk)

n=length(x(:,1))-1;
m=length(x(1,:))-1;



lsy=[dx:dx:dx*Ny]-dx*Ny/2-dx/2;
lsx=[dx:dx:dx*Nx]-dx*Nx/2-dx/2;
     
XL=max(lsx(Nx));
YL=max(lsy(Ny));

for i=1:kk:m
    hold off;
    figure (1)
    clf;
    hold on;
    plot(x(:,i),y(:,i),'.','Markersize',15);
    plot([-1.1*XL,1.2*XL],[1. *YL 1. *YL])
    plot([-1.1*XL,1.2*XL],[-1. *YL -1. *YL])
    plot([ 1.0*XL,1.0*XL],[-1.*YL  1.*YL],'r')
    
    %for j=1:n+1
    %    if x(j,i)>XL-dx
    %        plot(x(j,i)-2*XL+2*dx,y(j,i),'r*');
    %    end
    %end
    axis image;  
    axis([-1.11*XL 1.31*XL -1.01*YL 1.01*YL]);
    set(gca,'FontSize',30,'FontName','Times') 
    %set(gca,'XTick',[-2,0,2]);
    %set(gca,'yTick',[-2,0,2]);
    title(strcat('t=',num2str(t(i))),'FontSize',30,'FontName','Times');
    
    
 
    
    mov((i-1)/kk+1)=getframe(gcf);

end


end