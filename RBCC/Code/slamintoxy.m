function [forcx,forcy]=slamintoxy(m,sl,theta,dkap,slam)

slams=fd1(slam,m);



forcx=zeros(1,m);
forcy=zeros(1,m);


forcx(1,1:m)=(slam(1,1:m).*dkap(1,1:m).*sin(theta(1,1:m))...
-slams(1,1:m).*cos(theta(1,1:m))/sl);

forcy(1,1:m)=(-slam(1,1:m).*dkap(1,1:m).*cos(theta(1,1:m))...
-slams(1,1:m).*sin(theta(1,1:m))/sl);
end
