function [x,y]=vesiclereassign(m,x,y,ks)

    if ks==1
        x=x;
        y=y;
        
    elseif ks==-1
        
        
        n=round(length(x)/m);
        
        for i=1:n
            
            
            xc=mean(x((i-1)*m+1:i*m));
            x((i-1)*m+1:i*m)=-(x((i-1)*m+1:i*m)-xc)+xc;
            
            
            
            
            
            
            
            
        end
        
    end




end