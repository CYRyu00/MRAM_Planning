function [AM_com, AM_mass, AM_inertia] = get_inertia_ver2(shape ,m0, I0, d)
    [core_row, core_col] = find(shape == 2);
    [AMs_rows, AMs_cols] = find(shape ~= 0);
    
    m=length(AMs_cols);
    AM_com = [0;0;0]; AM_mass = m0*m; AM_inertia = zeros(3,3);
    for i =1:length(AMs_cols)
        r = [ (core_col - AMs_cols(i)) *-d ; (core_row - AMs_rows(i)) *d ;0];% p_j,core
        AM_com = AM_com - r*m0/AM_mass;
    end
    for i =1:length(AMs_cols)
        r = [ (core_col - AMs_cols(i)) *-d ; (core_row - AMs_rows(i)) *d ;0] ;%p_j to com
        AM_inertia = AM_inertia + I0 +m0*(r'*r*eye(3,3)-r*r');
    end
end

