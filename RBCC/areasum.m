function sum2=areasum(m,x,y)

    xcc=mean(x(1,1:m));
    ycc=mean(y(1,1:m));
    x(1:m)=x(1:m)-xcc;
    y(1:m)=y(1:m)-ycc;
    
    y(m+1)=y(1);
    x(m+1)=x(1);
 
    
    sum2=sum(abs(x(1,2:m+1).*( -y(1,1:m))+0.*(y(1,1:m)-y(1,2:m+1))...
    +x(1,1:m).*(y(1,2:m+1) )))/2;
    

end
