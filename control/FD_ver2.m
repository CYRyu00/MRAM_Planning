function qdd = FD_ver2(n, DH, mass, inertia, r_i_ci, g, q, qd, tau, F_ext)
    import casadi.*
     
    q = MX(q);
    qd = MX(qd);
    tau = MX(tau);
    F_ext = MX(F_ext);
    qdd_ = MX.zeros(n,1);

    h = newton_euler_inverse_dynamics(n, DH, mass, inertia, r_i_ci, g, q, qd, qdd_, -F_ext);

    M_matrix = MX.zeros(n, n);
    for i = 1:n
        g_ = MX([0; 0; 0]);
        qd_ = MX.zeros(n,1);
        F_ext_ = MX.zeros(6,1);
        qdd_ = MX.zeros(n,1);
        qdd_(i) = 1;

        M_matrix(i,:) = newton_euler_inverse_dynamics(n, DH, mass, inertia, r_i_ci, g_, q, qd_, qdd_, -F_ext_)';
    end

    qdd = M_matrix \ (tau - h);
end
