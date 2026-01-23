function [qdd, f_int, tau_int] = FD_cart_ver2(p_1, phi_1, pd_1, w_1, shape_pitch, r_g, U, M, k_obj, b_obj, l1, l2, n_o, n_am)
    e_1 = [1; 0; 0]; %e_2 = [0; 1; 0]; e_3 = [0; 0; 1];
    % kinematics
    qd = zeros(n_o + n_am*6, 1);
    p_o = p_1 - r_g + l1 * Ry(phi_1) * e_1;
    qd(1:3) = pd_1 + l1 * Ry(phi_1) * S(w_1) * e_1;
    
    phi(1) = phi_1;
    qd(4:6) = pd_1;
    pd_prev = pd_1;
    qd(n_o + (n_am + 1)*3 - 2 : n_o + (n_am + 1)*3) = w_1;
    
    for j = 2:n_am
        phi(j) = phi_1 + shape_pitch(j);
        qd(n_o + (n_am + j)*3 - 2 : n_o + (n_am + j)*3) = w_1;
        qd(1 + 3*j : 3 + 3*j) = pd_prev - l2*Ry(phi(j-1))*S(w_1)*e_1 - l1*Ry(phi(j))*S(w_1)*e_1 ;
    end
    
    for j = 2:n_am
        phi(j) = phi_1 + shape_pitch(j);
        qd(n_o + (n_am + j)*3 - 2 : n_o + (n_am + j)*3) = w_1;
        qd(1 + 3*j : 3 + 3*j) = pd_prev - l2*Ry(phi(j-1))*S(w_1)*e_1 - l1*Ry(phi(j))*S(w_1)*e_1 ;
    end

    % pfaffin constraints
    A_am = zeros((n_am-1)*6,n_am*6);
    Ad_am = zeros((n_am-1)*6,n_am*6);
    
    for j = 1:n_am - 1
        R1 = Ry(phi(j));
        R2 = Ry(phi(j+1));
        w1 = qd(n_o + (n_am + j)*3 - 2 : n_o + (n_am + j)*3);
        w2 = qd(n_o + (n_am + j)*3 + 1 : n_o + (n_am + j)*3 + 3);

        A_am(3*j-2: 3*j, 3*j-2: 3*j+3) = [eye(3), -eye(3)];
        A_am(3*j-2: 3*j, 3*(j+n_am)-2: 3*(j+n_am)+3) = [l2 * R1 * S(e_1), l1 * R2 * S(e_1)];
        A_am(3*(j+n_am-1)-2: 3*(j+n_am-1), 3*(j+n_am)-2: 3*(j+n_am)+3) = [eye(3), -eye(3)];
    
        Ad_am(3*j-2: 3*j, 3*(j+n_am)-2: 3*(j+n_am)+3) = [l2 * R1 * S(w1) * S(e_1), l1 * R2 * S(w2) * S(e_1)];
    end
    R1 = Ry(phi(1));
    
    A_o = [0 1 0; 0 0 1];
    A_o_am = [eye(3), -eye(3), zeros(3, (n_am-1)*3), l1 * R1 * S(e_1), zeros(3, (n_am-1)*3) ; ...
             zeros(3, (n_am + 1)*3), -eye(3), zeros(3, (n_am-1)*3)];
    A = [A_o, zeros(2, n_am * 6); A_o_am; zeros((n_am-1)*6, n_o), A_am];
    
    Ad_o = [0 0 0; 0 0 0];
    Ad_o_am = [zeros(3, (n_am+1)*3), l1 * R1 * S(w_1) * S(e_1), zeros(3, (n_am-1)*3) ; ...
             zeros(3, (2*n_am + 1)*3)];
    Ad = [Ad_o, zeros(2, n_am * 6); Ad_o_am; zeros((n_am-1)*6, n_o), Ad_am];
    
    A_dagger = ((A* (M\ A')) \ A) * inv(M);
    
    Cqd = zeros(6 * n_am + n_o) * qd;

    % Object input
    U_obj = - k_obj * p_o - b_obj * qd(1:3);
   
    F_int = A_dagger * ([U_obj; U] - Cqd) + ((A* (M\ A'))\ Ad) * qd;
    f_int = [F_int(3:5); F_int(9: n_am*3 +5) ]; % f 1 ~ n
    tau_int = [F_int(6:8); F_int(n_am*3 +6:end) ]; % tau 1 ~ n

    qdd = M\ ([U_obj ; U] -  A'* F_int);
end

function out = Ry(q)
    out = [cos(q), 0, sin(q);
           0, 1, 0;
           -sin(q), 0, cos(q)];
end