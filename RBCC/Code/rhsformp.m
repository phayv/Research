function rhs=rhsformp(m,Nx,Ny,dx,noos,rforcx,rforcy,pressure)



nk=Ny;
n=Nx;

rhs=zeros(1,3*Nx*Ny);


   nx= noos(1,:);
   ny= noos(2,:);
   a2= noos(3,:);
   a1= noos(4,:);
   d2= noos(5,:);
   d1= noos(6,:);

   rx1=dx+a1;
   rx2=a1;
   rx3=a2;
   rx4=dx+a2;
   
   ry1=dx+d1;
   ry2=d1;
   ry3=d2;
   ry4=dx+d2;
   
   rx1=rx1/dx;
   rx2=rx2/dx;   
   rx3=rx3/dx;  
   rx4=rx4/dx;
   
   ry1=ry1/dx;
   ry2=ry2/dx;   
   ry3=ry3/dx;  
   ry4=ry4/dx;   
   
   
   
   drx1=(5-2*rx1-sqrt(-7+12*rx1-4*rx1.^2))/8;
   drx2=(3-2*rx2+sqrt(1+4*rx2-4*rx2.^2))/8;
   drx3=(3-2*rx3+sqrt(1+4*rx3-4*rx3.^2))/8;   
   drx4=(5-2*rx4-sqrt(-7+12*rx4-4*rx4.^2))/8;

   dry1=(5-2*ry1-sqrt(-7+12*ry1-4*ry1.^2))/8;
   dry2=(3-2*ry2+sqrt(1+4*ry2-4*ry2.^2))/8;
   dry3=(3-2*ry3+sqrt(1+4*ry3-4*ry3.^2))/8;   
   dry4=(5-2*ry4-sqrt(-7+12*ry4-4*ry4.^2))/8;
   
   
 
for i=1:m
    if nx(i)<n-1 && ny(i)<nk-1 && nx(i)>=2 && ny(i)>=2
        
   
        rhs((nx(i)-2)*nk+ny(i)-1)=rhs((nx(i)-2)*nk+ny(i)-1)+drx1(i)*dry1(i)*rforcx(i);
        rhs((nx(i)-2)*nk+ny(i))  =rhs((nx(i)-2)*nk+ny(i))  +drx1(i)*dry2(i)*rforcx(i);
        rhs((nx(i)-2)*nk+ny(i)+1)=rhs((nx(i)-2)*nk+ny(i)+1)+drx1(i)*dry3(i)*rforcx(i);
        rhs((nx(i)-2)*nk+ny(i)+2)=rhs((nx(i)-2)*nk+ny(i)+2)+drx1(i)*dry4(i)*rforcx(i);

        rhs((nx(i)-1)*nk+ny(i)-1)=rhs((nx(i)-1)*nk+ny(i)-1)+drx2(i)*dry1(i)*rforcx(i);
        rhs((nx(i)-1)*nk+ny(i))  =rhs((nx(i)-1)*nk+ny(i))  +drx2(i)*dry2(i)*rforcx(i);
        rhs((nx(i)-1)*nk+ny(i)+1)=rhs((nx(i)-1)*nk+ny(i)+1)+drx2(i)*dry3(i)*rforcx(i);
        rhs((nx(i)-1)*nk+ny(i)+2)=rhs((nx(i)-1)*nk+ny(i)+2)+drx2(i)*dry4(i)*rforcx(i);    

        rhs((nx(i))*nk+ny(i)-1)=rhs((nx(i))*nk+ny(i)-1)+drx3(i)*dry1(i)*rforcx(i);
        rhs((nx(i))*nk+ny(i))  =rhs((nx(i))*nk+ny(i))  +drx3(i)*dry2(i)*rforcx(i);
        rhs((nx(i))*nk+ny(i)+1)=rhs((nx(i))*nk+ny(i)+1)+drx3(i)*dry3(i)*rforcx(i);
        rhs((nx(i))*nk+ny(i)+2)=rhs((nx(i))*nk+ny(i)+2)+drx3(i)*dry4(i)*rforcx(i);  

        rhs((nx(i)+1)*nk+ny(i)-1)=rhs((nx(i)+1)*nk+ny(i)-1)+drx4(i)*dry1(i)*rforcx(i);
        rhs((nx(i)+1)*nk+ny(i))  =rhs((nx(i)+1)*nk+ny(i))  +drx4(i)*dry2(i)*rforcx(i);
        rhs((nx(i)+1)*nk+ny(i)+1)=rhs((nx(i)+1)*nk+ny(i)+1)+drx4(i)*dry3(i)*rforcx(i);
        rhs((nx(i)+1)*nk+ny(i)+2)=rhs((nx(i)+1)*nk+ny(i)+2)+drx4(i)*dry4(i)*rforcx(i);      



        rhs((nx(i)-2)*nk+ny(i)-1+nk*n)=rhs((nx(i)-2)*nk+ny(i)-1+nk*n)+drx1(i)*dry1(i)*rforcy(i);
        rhs((nx(i)-2)*nk+ny(i)  +nk*n)=rhs((nx(i)-2)*nk+ny(i)  +nk*n)+drx1(i)*dry2(i)*rforcy(i);
        rhs((nx(i)-2)*nk+ny(i)+1+nk*n)=rhs((nx(i)-2)*nk+ny(i)+1+nk*n)+drx1(i)*dry3(i)*rforcy(i);
        rhs((nx(i)-2)*nk+ny(i)+2+nk*n)=rhs((nx(i)-2)*nk+ny(i)+2+nk*n)+drx1(i)*dry4(i)*rforcy(i);

        rhs((nx(i)-1)*nk+ny(i)-1+nk*n)=rhs((nx(i)-1)*nk+ny(i)-1+nk*n)+drx2(i)*dry1(i)*rforcy(i);
        rhs((nx(i)-1)*nk+ny(i)  +nk*n)=rhs((nx(i)-1)*nk+ny(i)  +nk*n)+drx2(i)*dry2(i)*rforcy(i);
        rhs((nx(i)-1)*nk+ny(i)+1+nk*n)=rhs((nx(i)-1)*nk+ny(i)+1+nk*n)+drx2(i)*dry3(i)*rforcy(i);
        rhs((nx(i)-1)*nk+ny(i)+2+nk*n)=rhs((nx(i)-1)*nk+ny(i)+2+nk*n)+drx2(i)*dry4(i)*rforcy(i);    

        rhs((nx(i))*nk+ny(i)-1+nk*n)=rhs((nx(i))*nk+ny(i)-1+nk*n)+drx3(i)*dry1(i)*rforcy(i);
        rhs((nx(i))*nk+ny(i)  +nk*n)=rhs((nx(i))*nk+ny(i)  +nk*n)+drx3(i)*dry2(i)*rforcy(i);
        rhs((nx(i))*nk+ny(i)+1+nk*n)=rhs((nx(i))*nk+ny(i)+1+nk*n)+drx3(i)*dry3(i)*rforcy(i);
        rhs((nx(i))*nk+ny(i)+2+nk*n)=rhs((nx(i))*nk+ny(i)+2+nk*n)+drx3(i)*dry4(i)*rforcy(i);  

        rhs((nx(i)+1)*nk+ny(i)-1+nk*n)=rhs((nx(i)+1)*nk+ny(i)-1+nk*n)+drx4(i)*dry1(i)*rforcy(i);
        rhs((nx(i)+1)*nk+ny(i)  +nk*n)=rhs((nx(i)+1)*nk+ny(i)  +nk*n)+drx4(i)*dry2(i)*rforcy(i);
        rhs((nx(i)+1)*nk+ny(i)+1+nk*n)=rhs((nx(i)+1)*nk+ny(i)+1+nk*n)+drx4(i)*dry3(i)*rforcy(i);
        rhs((nx(i)+1)*nk+ny(i)+2+nk*n)=rhs((nx(i)+1)*nk+ny(i)+2+nk*n)+drx4(i)*dry4(i)*rforcy(i);

        
    elseif ny(i)==nk-1 
        rhs((nx(i)-2)*nk+ny(i)-1)=rhs((nx(i)-2)*nk+ny(i)-1)+drx1(i)*dry1(i)*rforcx(i);
        rhs((nx(i)-2)*nk+ny(i))  =rhs((nx(i)-2)*nk+ny(i))  +drx1(i)*rforcx(i)/2;

        rhs((nx(i)-1)*nk+ny(i)-1)=rhs((nx(i)-1)*nk+ny(i)-1)+drx2(i)*dry1(i)*rforcx(i);
        rhs((nx(i)-1)*nk+ny(i))  =rhs((nx(i)-1)*nk+ny(i))  +drx2(i)*rforcx(i)/2;

        rhs((nx(i))*nk+ny(i)-1)=rhs((nx(i))*nk+ny(i)-1)+drx3(i)*dry1(i)*rforcx(i);
        rhs((nx(i))*nk+ny(i))  =rhs((nx(i))*nk+ny(i))  +drx3(i)*rforcx(i)/2;

        rhs((nx(i)+1)*nk+ny(i)-1)=rhs((nx(i)+1)*nk+ny(i)-1)+drx4(i)*dry1(i)*rforcx(i);
        rhs((nx(i)+1)*nk+ny(i))  =rhs((nx(i)+1)*nk+ny(i))  +drx4(i)*rforcx(i)/2;

        
        
        rhs((nx(i)-2)*nk+ny(i)-1+nk*n)=rhs((nx(i)-2)*nk+ny(i)-1+nk*n)+drx1(i)*dry1(i)*rforcy(i);
        rhs((nx(i)-2)*nk+ny(i)  +nk*n)=rhs((nx(i)-2)*nk+ny(i)  +nk*n)+drx1(i)/2*rforcy(i);

        rhs((nx(i)-1)*nk+ny(i)-1+nk*n)=rhs((nx(i)-1)*nk+ny(i)-1+nk*n)+drx2(i)*dry1(i)*rforcy(i);
        rhs((nx(i)-1)*nk+ny(i)  +nk*n)=rhs((nx(i)-1)*nk+ny(i)  +nk*n)+drx2(i)*rforcy(i)/2;

        rhs((nx(i))*nk+ny(i)-1+nk*n)=rhs((nx(i))*nk+ny(i)-1+nk*n)+drx3(i)*dry1(i)*rforcy(i);
        rhs((nx(i))*nk+ny(i)  +nk*n)=rhs((nx(i))*nk+ny(i)  +nk*n)+drx3(i)*rforcy(i)/2;

        rhs((nx(i)+1)*nk+ny(i)-1+nk*n)=rhs((nx(i)+1)*nk+ny(i)-1+nk*n)+drx4(i)*dry1(i)*rforcy(i);
        rhs((nx(i)+1)*nk+ny(i)  +nk*n)=rhs((nx(i)+1)*nk+ny(i)  +nk*n)+drx4(i)*rforcy(i)/2;

        
    elseif ny(i)==1
        
        rhs((nx(i)-2)*nk+ny(i)+1)=rhs((nx(i)-2)*nk+ny(i)+1)+drx1(i)/2*rforcx(i);
        rhs((nx(i)-2)*nk+ny(i)+2)=rhs((nx(i)-2)*nk+ny(i)+2)+drx1(i)*dry4(i)*rforcx(i);

        rhs((nx(i)-1)*nk+ny(i)+1)=rhs((nx(i)-1)*nk+ny(i)+1)+drx2(i)/2*rforcx(i);
        rhs((nx(i)-1)*nk+ny(i)+2)=rhs((nx(i)-1)*nk+ny(i)+2)+drx2(i)*dry4(i)*rforcx(i);    

        rhs((nx(i))*nk+ny(i)+1)=rhs((nx(i))*nk+ny(i)+1)+drx3(i)/2*rforcx(i);
        rhs((nx(i))*nk+ny(i)+2)=rhs((nx(i))*nk+ny(i)+2)+drx3(i)*dry4(i)*rforcx(i);  

        rhs((nx(i)+1)*nk+ny(i)+1)=rhs((nx(i)+1)*nk+ny(i)+1)+drx4(i)/2*rforcx(i);
        rhs((nx(i)+1)*nk+ny(i)+2)=rhs((nx(i)+1)*nk+ny(i)+2)+drx4(i)*dry4(i)*rforcx(i);      


        
        rhs((nx(i)-2)*nk+ny(i)+1+nk*n)=rhs((nx(i)-2)*nk+ny(i)+1+nk*n)+drx1(i)/2*rforcy(i);
        rhs((nx(i)-2)*nk+ny(i)+2+nk*n)=rhs((nx(i)-2)*nk+ny(i)+2+nk*n)+drx1(i)*dry4(i)*rforcy(i);
        
        rhs((nx(i)-1)*nk+ny(i)+1+nk*n)=rhs((nx(i)-1)*nk+ny(i)+1+nk*n)+drx2(i)/2*rforcy(i);
        rhs((nx(i)-1)*nk+ny(i)+2+nk*n)=rhs((nx(i)-1)*nk+ny(i)+2+nk*n)+drx2(i)*dry4(i)*rforcy(i);    

        rhs((nx(i))*nk+ny(i)+1+nk*n)=rhs((nx(i))*nk+ny(i)+1+nk*n)+drx3(i)/2*rforcy(i);
        rhs((nx(i))*nk+ny(i)+2+nk*n)=rhs((nx(i))*nk+ny(i)+2+nk*n)+drx3(i)*dry4(i)*rforcy(i);  

        rhs((nx(i)+1)*nk+ny(i)+1+nk*n)=rhs((nx(i)+1)*nk+ny(i)+1+nk*n)+drx4(i)/2*rforcy(i);
        rhs((nx(i)+1)*nk+ny(i)+2+nk*n)=rhs((nx(i)+1)*nk+ny(i)+2+nk*n)+drx4(i)*dry4(i)*rforcy(i);    
    end
end






for i=2:nk-1
    rhs(i+2*n*nk)=1/dx*pressure;
end



end



%         rhs((nx(i)-1)*nk+ny(i))=rhs((nx(i)-1)*nk+ny(i))+rx3(i)*ry3(i)*rforcy(i);
%         rhs((nx(i)-1)*nk+ny(i)+1)=rhs((nx(i)-1)*nk+ny(i)+1)+rx3(i)*ry2(i)*rforcy(i);
%         rhs(nx(i)*nk+ny(i))=rhs(nx(i)*nk+ny(i))+rx2(i)*ry3(i)*rforcy(i);
%         rhs(nx(i)*nk+ny(i)+1)=rhs(nx(i)*nk+ny(i)+1)+rx2(i)*ry2(i)*rforcy(i);
% 
%         rhs((nx(i)-1)*nk+ny(i)+nk*n)=rhs((nx(i)-1)*nk+ny(i)+nk*n)+rx3(i)*ry3(i)*rforcy(i);
%         rhs((nx(i)-1)*nk+ny(i)+1+nk*n)=rhs((nx(i)-1)*nk+ny(i)+1+nk*n)+rx3(i)*ry2(i)*rforcy(i);
%         rhs(nx(i)*nk+ny(i)+nk*n)=rhs(nx(i)*nk+ny(i)+nk*n)+rx2(i)*ry3(i)*rforcy(i);
%         rhs(nx(i)*nk+ny(i)+1+nk*n)=rhs(nx(i)*nk+ny(i)+1+nk*n)+rx2(i)*ry2(i)*rforcy(i);   
