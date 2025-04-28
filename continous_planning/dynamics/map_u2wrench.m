function wrench = map_u2wrench( u, rho,k,l, core , mu , r , d)
    import casadi.*
    wrench = MX.zeros(6,1);
    A = [r r -r -r;  -r r r -r;mu -mu mu -mu;0 0 0 0; 0 0 0 0; 1 1 1 1];

    for i=1:k
        for j=1:l
            r = [ (core(2) - j) *-d ; (core(1) - i) *d ;0];% p_j,core
            F_bi = rho(i,j)*A*u( 4*((i-1)*l + j)-3 : 4*((i-1)*l + j) );% [ moment; force]
            wrench(1:3) = wrench(1:3) + F_bi(1:3) - cross(r,F_bi(4:6));
            wrench(4:6) = wrench(4:6) + F_bi(4:6);
        end
    end
end
