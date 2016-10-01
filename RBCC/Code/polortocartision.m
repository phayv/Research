
function [rforcx,rforcy]=polortocartision(m,theta,forc,fforc)

rforcx=zeros(1,m);
rforcy=zeros(1,m);


rforcx(1,1:m) =forc(1,1:m).*sin(theta(1,1:m))+fforc(1,1:m).*cos(theta(1,1:m));
rforcy(1,1:m) =-forc(1,1:m).*cos(theta(1,1:m))+fforc(1,1:m).*sin(theta(1,1:m));

end
