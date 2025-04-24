function [AM_com, AM_mass, AM_inertia] = get_inertia_double(rho,k,l, core ,m0, I0, d)
    
    AM_com = zeros(3,1);
    AM_mass = 0; 
    AM_inertia = zeros(3,3);
    
    for i=1:k
        for j=1:l           
            AM_mass = AM_mass + rho(i,j)*m0;
        end
    end
    
    for i=1:k
        for j=1:l
            r = [ (core(2) - j) *-d ; (core(1) - i) *d ;0];% p_j,core
            AM_com = AM_com - rho(i,j)*r*m0/AM_mass;
        end
    end

    for i=1:k
        for j=1:l
            r = [ (core(2) - j) *-d ; (core(1) - i) *d ;0]+ AM_com; %p_j to com
            AM_inertia =  AM_inertia + rho(i,j)*( I0 +m0*(r'*r*eye(3,3)-r*r') );
        end
    end

end

