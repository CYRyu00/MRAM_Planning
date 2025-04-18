function wrench = map_u2wrench_double( u, shape , mu , r , d)
    
    wrench = zeros(6,1);
    [core_row, core_col] = find(shape == 2);
    [AMs_rows, AMs_cols] = find(shape ~= 0);
    
    A = [r r -r -r;  -r r r -r;mu -mu mu -mu;0 0 0 0; 0 0 0 0; 1 1 1 1];
    for i=1:length(AMs_rows)
        r = [ (core_col - AMs_cols(i)) *-d ; (core_row - AMs_rows(i)) *d ;0];% p_j,core
        
        F_bi = A*u(4*i-3:4*i);% [ moment; force]
        wrench(1:3) = wrench(1:3) + F_bi(1:3) - cross(r,F_bi(4:6));
        wrench(4:6) = wrench(4:6) + F_bi(4:6);
    end
end
