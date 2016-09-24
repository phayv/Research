function [noos,fx]=direct(m,nk,x,y)

fx=1/(nk-1);

%fist we need to compute the latic that works with the grid point
    lsy=[fx:fx:fx*nk]-fx*nk/2-fx/2;
    lsx=[fx:fx:fx*nk*4]-fx*nk*2-fx/2;
    a1=zeros(1,m);
    a2=zeros(1,m);
    d1=zeros(1,m);
    d2=zeros(1,m);
    nx=zeros(1,m);
    ny=zeros(1,m);

    %form the right hand side by the latice
    for i=1:m
        for jj=1:nk*4-1
            if lsx(jj)<= x(1,i) && lsx(jj+1)>x(1,i)
                nx(i)=jj;
                a1(i)=abs(lsx(nx(i)+1)-x(i));
                a2(i)=abs(x(i)-lsx(nx(i)));
            end
        end
        for jj=1:nk-1
            if lsy(jj)<= y(1,i) && lsy(jj+1)>y(1,i)
                ny(i)=jj;
                d1(i)=abs(lsy(ny(i)+1)-y(i));
                d2(i)=abs(y(i)-lsy(ny(i)));
            end
        end
    end

    noos=zeros(6,m);
    noos(1,:)=nx;
    noos(2,:)=ny;
    noos(3,:)=a1;
    noos(4,:)=a2;
    noos(5,:)=d1;
    noos(6,:)=d2;
    
    noos(:,m+1)=noos(:,1);
    
end
