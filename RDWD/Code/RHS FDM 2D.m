% RHS FDM 2D %
function F = RHS_FDM_2D_MorphogenGradient(N,h,f,u,v,x)

F = zeros(N^2,1);

for i = 1:N
    for j = 1:N;
        m = i + N*(j-1);
        F(m,1) = h^2*(f(u(m,1)) + v(x(i,j)));
    end
end

end
