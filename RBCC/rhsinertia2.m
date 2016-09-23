function rhs=rhsinertia2(m,L,vel,dt,rho)


n=m*L;

rhs=zeros(1,3*m*n);
vel2=zeros(2*m*n,1);

vel2(1:n*3/4*m)=vel(n/4*m+1:n*m);
vel2(n*m+1:n*m+n*3/4*m)=vel(n*m+n/4*m+1:2*n*m);


rhs(1,1:2*m*n)=-vel2(1:2*m*n)'*rho/dt;


    for i=1:n
        rhs(1,i*m)=0;
        rhs(1,i*m+m*n)=0;
        rhs(1,(i-1)*m+1)=0;
        rhs(1,(i-1)*m+1+m*n)=0;
    end   

    for i=2:m-1
        rhs(1,i)=0;
    end
    

end
