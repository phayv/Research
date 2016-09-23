function noos=directperid(m,Nx,Ny,x,y,dx)

lm=length(x);
noos=zeros(6,lm);

%fist we need to compute the latic that works with the grid point
    lsy=[dx:dx:dx*Ny]-dx*Ny/2-dx/2;
    lsx=[dx:dx:dx*Nx]-dx*Nx/2-dx/2;
    

    for i=1:lm
        if x(i)<lsx(Nx-2)
            for jj=1:Nx-1
                if lsx(jj)<= x(1,i) && lsx(jj+1)>x(1,i)
                    nx(i)=jj;
                    a1(i)=abs(lsx(nx(i)+1)-x(i));
                    a2(i)=abs(x(i)-lsx(nx(i)));
                end
            end
            for jj=1:Ny-1
                if lsy(jj)<= y(1,i) && lsy(jj+1)>y(1,i)
                    ny(i)=jj;
                    d1(i)=abs(lsy(ny(i)+1)-y(i));
                    d2(i)=abs(y(i)-lsy(ny(i)));
                end
            end
        else 
           x(i)=x(i)-lsx(Nx-2)*2;
            for jj=1:Nx-1
                if lsx(jj)<= x(1,i) && lsx(jj+1)>x(1,i)
                    nx(i)=jj;
                    a1(i)=abs(lsx(nx(i)+1)-x(i));
                    a2(i)=abs(x(i)-lsx(nx(i)));
                end
            end
            for jj=1:Ny-1
                if lsy(jj)<= y(1,i) && lsy(jj+1)>y(1,i)
                    ny(i)=jj;
                    d1(i)=abs(lsy(ny(i)+1)-y(i));
                    d2(i)=abs(y(i)-lsy(ny(i)));
                end
            end
        end
    end


        noos(1,:)=nx;
        noos(2,:)=ny;
        noos(3,:)=a1;
        noos(4,:)=a2;
        noos(5,:)=d1;
        noos(6,:)=d2;

end
