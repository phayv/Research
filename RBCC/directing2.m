function noos=directing2(m,x,y,dx)
    n=length(x)-1;

%fist we need to compute the latic that works with the grid point
    ls=[dx:dx:dx*m]-dx*m/2-dx/2;
    a1=zeros(1,n);
    a2=zeros(1,n);
    d1=zeros(1,n);
    d2=zeros(1,n);
    nx=zeros(1,n);
    ny=zeros(1,n);

    %form the right hand side by the latice
    for i=1:n
        for jj=1:m-1
            if ls(jj)<= x(1,i) && ls(jj+1)>x(1,i)
                nx(i)=jj;
                a1(i)=abs(ls(nx(i)+1)-x(i));
                a2(i)=abs(x(i)-ls(nx(i)));
            end
            if ls(jj)<= y(1,i) && ls(jj+1)>y(1,i)
                ny(i)=jj;
                d1(i)=abs(ls(ny(i)+1)-y(i));
                d2(i)=abs(y(i)-ls(ny(i)));
            end
        end
    end

    noos=zeros(6,n);
    noos(1,:)=nx;
    noos(2,:)=ny;
    noos(3,:)=a1;
    noos(4,:)=a2;
    noos(5,:)=d1;
    noos(6,:)=d2;
    
    noos(:,n+1)=noos(:,1);
    

end
