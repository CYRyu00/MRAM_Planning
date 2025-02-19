function tau = newton_euler_inverse_dynamics(n, DH, mass, inertia, r_i_ci, g, q, qd, qdd, F_ext)
    import casadi.*
    
    
    q = MX(q);
    qd = MX(qd);
    qdd = MX(qdd);
    g = MX(g);
    F_ext = MX(F_ext);

    
    z0 = MX([0; 0; 1]);
    f = cell(n+1, 1);
    mu = cell(n+1, 1);
    tau = MX.zeros(n,1);  
    
    w = cell(n+1, 1); 
    alpha = cell(n+1, 1); 
    a = cell(n+1, 1);     
    a_e = cell(n+1,1);
    q = [q; 0];

    % Forward recursion
    for i = 1:n
        R = R_DH(DH(i,:), q(i)); 
        p = p_DH(DH(i+1,:), q(i+1));

        bi = z0;
        if i == 1
            w{i} = bi * qd(i);
            alpha{i} = bi * qdd(i);
            a{i} = R' * (-g) + cross(alpha{i}, r_i_ci{i}) ...
                   + cross(w{i}, cross(w{i}, r_i_ci{i}));
            a_e{i} = R' * (-g) + cross(alpha{i}, p) + cross(w{i}, cross(w{i}, p));
        else
            w{i} = R' * w{i-1} + bi * qd(i);
            alpha{i} = R' * alpha{i-1} + bi * qdd(i) + R' * cross(w{i-1}, bi * qd(i));
            a{i} = R' * a_e{i-1} + cross(alpha{i}, r_i_ci{i}) ...
                   + cross(w{i}, cross(w{i}, r_i_ci{i}));
            a_e{i} = R' * a_e{i-1} + cross(alpha{i}, p) + cross(w{i}, cross(w{i}, p));
        end 
    end

    % Backward recursion
    f{n+1} = F_ext(4:6);
    mu{n+1} = F_ext(1:3);

    for i = n:-1:1
        R = R_DH(DH(i+1,:), q(i+1)); 
        p = p_DH(DH(i+1,:), q(i+1));
        r_i1_ci = R' * (r_i_ci{i} - p);

        f{i} = R * f{i+1} + mass{i} * a{i};
        mu{i} = R * mu{i+1} - cross(f{i}, r_i_ci{i}) + R * cross(f{i+1}, r_i1_ci) ...
                + inertia{i} * alpha{i} + cross(w{i}, inertia{i} * w{i});
        tau(i) = mu{i}' * z0;
    end
end

function R = R_DH(dh, q)
    import casadi.*
    c1 = cos(dh(1)); c2 = cos(dh(4) + q); 
    s1 = sin(dh(1)); s2 = sin(dh(4) + q);
    
    R = [c2    -s2     0 ; 
         c1*s2 c1*c2 -s1 ; 
         s1*s2 c2*s1  c1 ]; 
end

function p = p_DH(dh, q)
    import casadi.*
    p = MX([dh(2); -sin(dh(1)) * dh(3); cos(dh(1)) * dh(3)]);
end
