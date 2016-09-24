function [uxs,uys]=areaconserve(ux,uy,thetas,Nv)

Mn=length(ux)/Nv;
 
for i=1:Nv
 
     uxI=ux(1,(i-1)*Mn+1:i*Mn);
     uyI=uy(1,(i-1)*Mn+1:i*Mn);
     theta(1,:)=thetas(i,:);

     [un,utt]=cartisionintopolar(Mn,theta,uxI,uyI) ;
     un=un-mean(un);  
     [uxk,uyk]=polortocartision(Mn,theta,un,utt);
     uxs(1,(i-1)*Mn+1:i*Mn)=uxk;
     uys(1,(i-1)*Mn+1:i*Mn)=uyk;
    

end


end
