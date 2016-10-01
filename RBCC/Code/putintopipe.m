function [x0,y0]=putintopipe(x,y,Nv)
    m=length(x);
    if mod(Nv,2)==0
        x=x-1.3/2;
    end
    LR=-1;
    for i=1:Nv
        x0((i-1)*m+1:i*m)=x+LR^(i)*1.3*floor(i/2);
        y0((i-1)*m+1:i*m)=y;
    end
end
