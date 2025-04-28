%Forward Dynamics using Newton - Euler Inverse Dynamics
function qdd = FD(n, A, M, mass, inertia, F_tip, q, qd, tau, g)
    qdd_ = zeros(n,1);
    % h(qd,qdd) + J'*F_tip
    h = newton_euler_ID(n, A, M, mass, inertia, F_tip, q, qd, qdd_, g);
    
    M_matrix = zeros(n,n);
    for i=1:n
        g_ = [0;0;0];
        qd_ = zeros(n,1);
        F_tip_ = [0;0;0;0;0;0];
        qdd_ = zeros(n,1); qdd_(i) = 1;
        M_matrix(i,:) = newton_euler_ID(n, A, M, mass, inertia, F_tip_, q, qd_, qdd_, g_);
    end
    
    
    qdd = M_matrix\(tau - h);
    
end