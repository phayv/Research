
function [run,rut]=cartisionintopolar(m,theta,ux,uy)

run=zeros(1,m);
rut=zeros(1,m);


run(1,1:m)=ux(1,1:m).*sin(theta(1,1:m))-uy(1,1:m).*cos(theta(1,1:m));
rut(1,1:m)=ux(1,1:m).*cos(theta(1,1:m))+uy(1,1:m).*sin(theta(1,1:m));

end

