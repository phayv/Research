function [x,y,sl]=initiall(shortax,ngrid)

     tol=10^(-12);
     N=(0:1:ngrid);
     x0=cos(2*pi*N/ngrid);
     y0=shortax*sin(2*pi*N/ngrid);
     
     dx=fd1(x0,ngrid);
     dy=fd1(y0,ngrid);
     
     slo=sqrt(dx.^2+dy.^2);
     slo(ngrid+1)=slo(1);
     alpha=arc(slo,ngrid);

     x=cos(pi*2*alpha);
     y=shortax*sin(pi*2*alpha);
     
     sxn=fd1(x,ngrid);
     syn=fd1(y,ngrid);
     sl=sum(sqrt(sxn(1:ngrid).^2+syn(1:ngrid).^2));
     sl=sl/ngrid;
end
 