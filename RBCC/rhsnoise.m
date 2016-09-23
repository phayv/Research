function rhs=rhsnoise(m,Nx,Ny,dx,Noi)


n=Nx;
m=Ny;


rhs=zeros(1,3*m*n);

rhs(1,1:2*m*n)=random('Normal',0,Noi/dx,1,2*m*n);



for i=1:n
    rhs(1+(i-1)*m)=0;
    rhs(i*m)=0;
    rhs(1+(i-1)*m+m*n)=0;
    rhs(i*m+m*n)=0;

end

for i=2:m-1
    rhs(i)=0;
    rhs(m*n-i+1)=0;
    rhs(m*n+i)=0;
    rhs(2*m*n-i+1)=0;

end

end
