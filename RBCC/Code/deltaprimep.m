
function [ux,uy]=deltaprimep(Nx,Ny,dx,noos,unn,utt)

nk=Ny;
m=length(noos(1,:));
ux=zeros(1,m);
uy=zeros(1,m);

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

    if nx(i)<Nx-1 && ny(i)<nk-1 && nx(i)>=2 && ny(i)>=2
        
        ux(1,i)=( unn((nx(i)-2)*nk+ny(i)-1)*drx1(i)*dry1(i)+unn((nx(i)-2)*nk+ny(i))*drx1(i)*dry2(i)+...
        unn((nx(i)-2)*nk+ny(i)+1)*drx1(i)*dry3(i)+unn((nx(i)-2)*nk+ny(i)+2)*drx1(i)*dry4(i)+...
        unn((nx(i)-1)*nk+ny(i)-1)*drx2(i)*dry1(i)+unn((nx(i)-1)*nk+ny(i))*drx2(i)*dry2(i)+...
        unn((nx(i)-1)*nk+ny(i)+1)*drx2(i)*dry3(i)+unn((nx(i)-1)*nk+ny(i)+2)*drx2(i)*dry4(i)+...
        unn((nx(i))*nk+ny(i)-1)*drx3(i)*dry1(i)+unn((nx(i))*nk+ny(i))*drx3(i)*dry2(i)+...
        unn((nx(i))*nk+ny(i)+1)*drx3(i)*dry3(i)+unn((nx(i))*nk+ny(i)+2)*drx3(i)*dry4(i)+... 
        unn((nx(i)+1)*nk+ny(i)-1)*drx4(i)*dry1(i)+unn((nx(i)+1)*nk+ny(i))*drx4(i)*dry2(i)+...
        unn((nx(i)+1)*nk+ny(i)+1)*drx4(i)*dry3(i)+unn((nx(i)+1)*nk+ny(i)+2)*drx4(i)*dry4(i));
    
        uy(1,i)=( utt((nx(i)-2)*nk+ny(i)-1)*drx1(i)*dry1(i)+utt((nx(i)-2)*nk+ny(i))*drx1(i)*dry2(i)+...
        utt((nx(i)-2)*nk+ny(i)+1)*drx1(i)*dry3(i)+utt((nx(i)-2)*nk+ny(i)+2)*drx1(i)*dry4(i)+...
        utt((nx(i)-1)*nk+ny(i)-1)*drx2(i)*dry1(i)+utt((nx(i)-1)*nk+ny(i))*drx2(i)*dry2(i)+...
        utt((nx(i)-1)*nk+ny(i)+1)*drx2(i)*dry3(i)+utt((nx(i)-1)*nk+ny(i)+2)*drx2(i)*dry4(i)+...
        utt((nx(i))*nk+ny(i)-1)*drx3(i)*dry1(i)+utt((nx(i))*nk+ny(i))*drx3(i)*dry2(i)+...
        utt((nx(i))*nk+ny(i)+1)*drx3(i)*dry3(i)+utt((nx(i))*nk+ny(i)+2)*drx3(i)*dry4(i)+... 
        utt((nx(i)+1)*nk+ny(i)-1)*drx4(i)*dry1(i)+utt((nx(i)+1)*nk+ny(i))*drx4(i)*dry2(i)+...
        utt((nx(i)+1)*nk+ny(i)+1)*drx4(i)*dry3(i)+utt((nx(i)+1)*nk+ny(i)+2)*drx4(i)*dry4(i));
        
%     elseif ny(i)==nk-1
%         ux(1,i)=( unn((nx(i)-2)*nk+ny(i)-1)*drx1(i)*dry1(i)+unn((nx(i)-2)*nk+ny(i))*drx1(i)/2+...
%         unn((nx(i)-1)*nk+ny(i)-1)*drx2(i)*dry1(i)+unn((nx(i)-1)*nk+ny(i))*drx2(i)/2+...
%         unn((nx(i))*nk+ny(i)-1)*drx3(i)*dry1(i)+unn((nx(i))*nk+ny(i))*drx3(i)/2+...
%         unn((nx(i)+1)*nk+ny(i)-1)*drx4(i)*dry1(i)+unn((nx(i)+1)*nk+ny(i))*drx4(i)/2);
%     
%         uy(1,i)=( utt((nx(i)-2)*nk+ny(i)-1)*drx1(i)*dry1(i)+utt((nx(i)-2)*nk+ny(i))*drx1(i)/2+...
%         utt((nx(i)-1)*nk+ny(i)-1)*drx2(i)*dry1(i)+utt((nx(i)-1)*nk+ny(i))*drx2(i)/2+...
%         utt((nx(i))*nk+ny(i)-1)*drx3(i)*dry1(i)+utt((nx(i))*nk+ny(i))*drx3(i)/2+...
%         utt((nx(i)+1)*nk+ny(i)-1)*drx4(i)*dry1(i)+utt((nx(i)+1)*nk+ny(i))*drx4(i)/2);        
%     elseif ny(i)==1
%         ux(1,i)=(unn((nx(i)-2)*nk+ny(i)+1)*drx1(i)/2+unn((nx(i)-2)*nk+ny(i)+2)*drx1(i)*dry4(i)+...
%         unn((nx(i)-1)*nk+ny(i)+1)*drx2(i)/2+unn((nx(i)-1)*nk+ny(i)+2)*drx2(i)*dry4(i)+...
%         unn((nx(i))*nk+ny(i)+1)*drx3(i)/2+unn((nx(i))*nk+ny(i)+2)*drx3(i)*dry4(i)+... 
%         unn((nx(i)+1)*nk+ny(i)+1)*drx4(i)/2+unn((nx(i)+1)*nk+ny(i)+2)*drx4(i)*dry4(i));
%     
%         uy(1,i)=(utt((nx(i)-2)*nk+ny(i)+1)*drx1(i)/2+utt((nx(i)-2)*nk+ny(i)+2)*drx1(i)*dry4(i)+...
%         utt((nx(i)-1)*nk+ny(i)+1)*drx2(i)/2+utt((nx(i)-1)*nk+ny(i)+2)*drx2(i)*dry4(i)+...
%         utt((nx(i))*nk+ny(i)+1)*drx3(i)/2+utt((nx(i))*nk+ny(i)+2)*drx3(i)*dry4(i)+... 
%         utt((nx(i)+1)*nk+ny(i)+1)*drx4(i)/2+utt((nx(i)+1)*nk+ny(i)+2)*drx4(i)*dry4(i));     
    end
    
end

end


%     ux(1,i)=1/dx^2*((a2(i)*d2(i))*unn((nx(i)-1)*nk+ny(i))+(a2(i)*d1(i))*unn((nx(i)-1)*nk+ny(i)+1)...
%     +(a1(i)*d2(i))*unn(nx(i)*nk+ny(i))+(a1(i)*d1(i))*unn(nx(i)*nk+ny(i)+1));
%     uy(1,i)=1/dx^2*((a2(i)*d2(i))*utt((nx(i)-1)*nk+ny(i))+(a2(i)*d1(i))*utt((nx(i)-1)*nk+ny(i)+1)...
%     +(a1(i)*d2(i))*utt(nx(i)*nk+ny(i))+(a1(i)*d1(i))*utt(nx(i)*nk+ny(i)+1));
