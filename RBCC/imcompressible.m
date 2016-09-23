
function vel=imcompressible(m,dx,sl,theta,dkap,noos,ps)

n=m;

unn=ps(1,1:n^2);
utt=ps(1,n^2+1:2*n^2);

%put unn into ux uy

[ux,uy]=deltaprime(m,dx,noos,unn,utt);

%put cartistion into polar

[run,rut]=cartisionintopolar(m,theta,ux,uy);

ruts=fd1(rut,m);
vel=(ruts(1,1:m)/sl+run(1,1:m).*dkap(1,1:m));
vel(1,m)=sum(run);
end
