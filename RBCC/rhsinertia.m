function rhs=rhsinertia(m,L,vel,dt,rho)


n=m*L;

rhs=zeros(1,3*m*n);

rhs(1,1:2*m*n)=-vel(1:2*m*n)'*rho/dt;


    for i=1:n
        rhs(1,i*m)=0;
        rhs(1,i*m+m*n)=0;
        rhs(1,(i-1)*m+1)=0;
        rhs(1,(i-1)*m+1+m*n)=0;
    end   


end
