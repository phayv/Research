function x1=perioddet(Nv,x,Nx,dx)
  x1=x;  
  lm=length(x);
  m=round(lm/Nv);
  lsx=[dx:dx:dx*Nx]-dx*Nx/2-dx/2;
  for i=1:Nv
      xI=x1(1,(i-1)*m+1:i*m);
      if min(xI)>lsx(Nx-2)
          x1(1,(i-1)*m+1:i*m)=x1(1,(i-1)*m+1:i*m)-2*lsx(Nx-2);
      end
  end

end
