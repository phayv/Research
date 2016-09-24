% D RHS FDM 2D %
function Fp = derivativeRHS_FDM_2D_Morphogen(N,h,fp,u)

Fp = zeros(N^2,1);

for i = 1:N^2
    Fp(i,1) = h^2*fp(u(i,1));
end

end
