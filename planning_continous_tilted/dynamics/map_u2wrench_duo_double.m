function wrench = map_u2wrench_duo_double(u, rho, k, l, core, mu, r, d, theta)
    wrench = zeros(6,1);
    
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
    
    for i=1:k
        for j=1:l
            r = [ (core(2) - j) * -2 *d ; (core(1) - i) *d ;0];% p_j,core
            F_bi = rho(i,j) * A_duo * u( 8*((i-1)*l + j) - 7 : 8*((i-1)*l + j) );% [ moment; force]
            wrench(1:3) = wrench(1:3) + F_bi(1:3) - cross(r,F_bi(4:6));
            wrench(4:6) = wrench(4:6) + F_bi(4:6);
        end
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
