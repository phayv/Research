function W=formWp(Nx,Ny,dx,mu)

    L=Nx;
    n=Ny;
    M=Nx*Ny;

    coef1=mu/dx^2;
    coef2=1/2/dx;
    
    
    
    % u, v
    K1=[1:1:2*M];
    L1=[1:1:2*M];
    S1=-4*coef1*ones(1,2*M);
    K2=[2:1:2*M];
    L2=[1:1:2*M-1];
    S2=coef1*ones(1,2*M-1);
    K3=[1:1:2*M-1];
    L3=[2:1:2*M];
    S3=coef1*ones(1,2*M-1);
    K4=[1:1:2*M-n];
    L4=[n+1:1:2*M];
    S4=coef1*ones(1,2*M-n);
    K5=[n+1:1:2*M];
    L5=[1:1:2*M-n];
    S5=coef1*ones(1,2*M-n);
    % P
    K6=[1:1:M-n];
    L6=[2*M+n+1:3*M];
    S6=-coef2*ones(1,M-n);
    K7=[n+1:1:M];
    L7=[2*M+1:3*M-n];
    S7=coef2*ones(1,M-n);
    K8=[M+1:1:2*M-1];
    L8=[2*M+2:3*M];
    S8=-coef2*ones(1,M-1);
    K9=[M+2:1:2*M];
    L9=[2*M+1:3*M-1];
    S9=coef2*ones(1,M-1);
    % imcompressible
    W=sparse([L1,L2,L3,L4,L5,L6,L7,L8,L9,K6,K7,K8,K9],...
    [K1,K2,K3,K4,K5,K6,K7,K8,K9,L6,L7,L8,L9],...
    [S1,S2,S3,S4,S5,S6,S7,S8,S9,S6,S7,S8,S9]);
    %x direction
    for j=2:n-1
        % the first column coef1*(-2*v(j)+v(j+1)+v(j-1))...
        %    +coef1*(-v(j)+v(j+n))...
        %    -coef2*(-3*p(j)+4*p(j+n)-p(j+2*n));
        W(j,j)=-3*coef1;
        W(j,2*M+j)=3*coef2;
        W(j,2*M+j+n)=-4*coef2;
        W(j,2*M+j+2*n)=coef2;
        % the last column%
        % vel(j+M-n)=coef1*(-2*v(j+M-n)+v(j+1+M-n)+v(j-1+M-n))...
        %    +coef1*(-v(j+M-n)+v(j+M-2*n))...
        %    +coef2*(-3*p(j+M-n)+4*p(j+M-2*n)-p(j+M-3*n));
        W(j+M-n,j+M-n)=-3*coef1;
        W(j+M-n,j+M)=0;
        W(j+M-n,3*M+j-n)=-3*coef2;
        W(j+M-n,3*M+j-2*n)=4*coef2;
        W(j+M-n,3*M+j-3*n)=-coef2;
    end
    for j=2:L
        % the bot row%
        %vel(1+(j-1)*n)=-4*coef1*v(1+(j-1)*n);
        W(1+(j-1)*n,1+(j-1)*n+1)=0;
        W(1+(j-1)*n,1+(j-1)*n+n)=0;
        W(1+(j-1)*n,1+(j-1)*n-1)=0;
        W(1+(j-1)*n,1+(j-1)*n-n)=0;
        W(1+(j-2)*n,1+(j-2)*n+2*M+n)=0;
        W(1+(j-1)*n,1+(j-1)*n+2*M-n)=0;
        % the top row%
        %vel(j*n)=-4*coef1*v(j*n);
        W((j-1)*n,(j-1)*n+1)=0;
        W((j-1)*n,(j-1)*n+n)=0;
        W((j-1)*n,(j-1)*n-1)=0;
        W(j*n,j*n-n)=0;
        W((j-1)*n,(j-1)*n+2*M+n)=0;
        W((j-1)*n,(j-1)*n+2*M-n)=0;
    end
    W(j*n,(j-1)*n+2*M)=0;
    W(1,2)=0;
    W(1,n+1)=0;
    W(1,2*M+1+n)=0;
    W(M,M+1)=0;
    W(M,M-1)=0;
    W(M,M+n)=0;
    W(M,M-n)=0;
    W(M,3*M-1)=0;
    %%%y direction
    for j=2:n-1
        %the first column
        %             vel(j+M)=coef1*(-2*u(j)+u(j+1)+u(j-1))...
        %                 +coef1*(-u(j)+u(j+n))...
        %                -coef2*(p(j+1)-p(j-1));
        W(j+M,j+M)=-3*coef1;
        W(j+M,j+M-n)=0;
        %the last column
        %             vel(j+2*M-n)=coef1*(-2*u(j+M-n)+u(j+1+M-n)+u(j-1+M-n))...
        %                 +coef1*(-u(j+M-n)+u(j+M-2*n))...
        %                -coef2*(p(j+M-n+1)-p(j+M-n-1));
        W(j+2*M-n,j+2*M-n)=-3*coef1;
        % W(j+2*M-n,j+2*M-2*n)=coef1;
        % W(j+2*M-n,j+2*M-2*n)=0;
    end
    for j=1:L
        % the bot row
        %vel(M+(j-1)*n+1)=-4*coef1*u((j-1)*n+1);
        W(M+(j-1)*n+1,M+(j-1)*n+2)=0;
        W(M+(j-1)*n+1,M+(j-1)*n)=0;
        W(M+(j-1)*n+1,M+(j-1)*n+1-n)=0;
        W(M+(j-2)*n+1,M+(j-2)*n+1+n)=0;
        % the top row
        %vel(M+j*n)=-4*coef1*u(j*n);
        W(M+j*n,M+j*n+1)=0;
        W(M+j*n,M+j*n-1)=0;
        W(M+j*n,M+j*n+n)=0;
        W(M+j*n,M+j*n-n)=0;
    end
    for j=1:L
        W((j-1)*n+M+1,(j-1)*n+2*M+2)=0;
        W((j-1)*n+M+1,(j-1)*n+2*M)=0;
        W(j*n+M,j*n+2*M-1)=0;
        W((j-1)*n+M,(j-1)*n+2*M+1)=0;
    end
    for j=1:L
        %%the bot row
        %            vel((j-1)*n+1+2*M)=coef2*u((j-1)*n+2);
        W((j-1)*n+1+2*M,(j-1)*n+2+M)=coef2;
        W((j-1)*n+1+2*M,j*n+1)=0;
        W((j-1)*n+1+2*M,M+(j-1)*n)=0;
        %%the top row
        %            vel(j*n+2*M)=coef2*u(j*n-1);
        W(j*n+2*M,j*n-1+M)=coef2;
        W(j*n+2*M,j*n+n)=0;
        W(j*n+2*M,j*n+M+1)=0;
    end
    for j=2:L
        W((j-1)*n+1+2*M,(j-2)*n+1)=0;
        W(j*n+2*M,j*n-n)=0;
    end
    for j=2:n-1
    %the first column
        W(j+2*M,j+2*M)=coef2;
        W(j+2*M,j+n)=0;
        W(j+3*M-n,j+3*M-n)=coef2;
        W(j+3*M-n,j+M-2*n)=0;
        W(j+2*M,j+M-1)=0;
        W(j+2*M,j+M+1)=0;
        W(j+3*M-n,j+2*M-n-1)=0;
        W(j+3*M-n,j+2*M-n+1)=0;
    end

    
    
    
    

end
