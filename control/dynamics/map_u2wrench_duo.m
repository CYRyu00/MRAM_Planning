function wrench = map_u2wrench_duo(u, shape, mu, r, d, theta)
    import casadi.*
    wrench = MX.zeros(6,1);
    [core_row, core_col] = find(shape == 2);
    [AMs_rows, AMs_cols] = find(shape ~= 0);
    
    s = sin(theta); c = cos(theta);
    A_theta = [mu*s/sqrt(2) - r*c, mu*s/sqrt(2) - r*c, - mu*s/sqrt(2) + r*c, - mu*s/sqrt(2) + r*c;
         - mu*s/sqrt(2) + r*c, mu*s/sqrt(2) - r*c, mu*s/sqrt(2) - r*c, - mu*s/sqrt(2) + r*c;
         mu*c + sqrt(2)*r*s, - mu*c - sqrt(2)*r*s, mu*c + sqrt(2)*r*s, - mu*c - sqrt(2)*r*s;
         s/sqrt(2), - s/sqrt(2), - s/sqrt(2), s/sqrt(2);
         - s/sqrt(2), - s/sqrt(2), s/sqrt(2), s/sqrt(2);
         c, c, c, c];
    
    p1 = [-0.5 * d; 0; 0]; p2 = [0.5 * d; 0; 0]; % p(:,3) = [D; 0; 0];
    R1 = eye(3,3); R2 = [0 -1 0;1 0 0; 0 0 1]; % R{3} = eye(3,3);
    
    A_duo = [Ad(R1, p1) * A_theta, Ad(R2, p2) * A_theta];

    for i=1:length(AMs_rows)
        r = [ (core_col - AMs_cols(i)) * -2 *d ; (core_row - AMs_rows(i)) *d ;0];% p_j,core
        
        F_bi = A_duo * u(8*i -7 : 8*i);% [moment; force]
        wrench(1:3) = wrench(1:3) + F_bi(1:3) - cross(r,F_bi(4:6));
        wrench(4:6) = wrench(4:6) + F_bi(4:6);
    end
end
function out = Ad(R, p)
out = [R zeros(3,3); S(p)*R R];
end
function out = S(p)
out = [0 -p(3) p(2);
       p(3) 0 -p(2);
       -p(2) p(1) 0];
end

