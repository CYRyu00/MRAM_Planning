function [qdd, f_in, tau_in] = FD_cart(p_o, phi, qd, U, M, k_obj, b_obj, l1, l2, n_o, n_am)
    e_1 = [1; 0; 0]; %e_2 = [0; 1; 0]; e_3 = [0; 0; 1];
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
    w1 = qd(n_o + (n_am + 1)*3 - 2 : n_o + (n_am + 1)*3);
    
    A_o = [0 1 0; 0 0 1];
    A_o_am = [eye(3), -eye(3), zeros(3, (n_am-1)*3), l1 * R1 * S(e_1), zeros(3, (n_am-1)*3) ; ...
             zeros(3, (n_am + 1)*3), eye(3), zeros(3, (n_am-1)*3)];
    A = [A_o, zeros(2, n_am * 6); A_o_am; zeros((n_am-1)*6, n_o), A_am];
    
    Ad_o = [0 0 0; 0 0 0];
    Ad_o_am = [zeros(3, (n_am+1)*3), l1 * R1 * S(w1) * S(e_1), zeros(3, (n_am-1)*3) ; ...
             zeros(3, (2*n_am + 1)*3)];
    Ad = [Ad_o, zeros(2, n_am * 6); Ad_o_am; zeros((n_am-1)*6, n_o), Ad_am];
    
    A_dagger = A' * inv((A* (M\ A')));
    
    A_f_t = zeros(n_am * 3, n_am * 3);
    A_f_r = zeros(n_am * 3, n_am * 3);
    A_tau_r = zeros(n_am * 3, n_am * 3);
    for j = 1 : n_am % f2, f3, f4 ...
        Rj = Ry(phi(j));
        if j == n_am
            A_f_t(3*j-2:3*j, 3*j-2:3*j) = eye(3);
            A_f_r(3*j-2:3*j, 3*j-2:3*j) = l1*S(Rj*e_1);
            A_tau_r(3*j-2:3*j, 3*j-2:3*j) = eye(3);
        else
            A_f_t(3*j-2:3*j, 3*j-2:3*j+3) = [eye(3), -eye(3)];
            A_f_r(3*j-2:3*j, 3*j-2:3*j+3) = [l1*S(Rj*e_1), -l2*S(Rj*e_1)];
            A_tau_r(3*j-2:3*j, 3*j-2:3*j+3) = [eye(3), -eye(3)];
        end
    end

    % Object input
    U_obj = - k_obj * p_o - b_obj * qd(1:3);
   
    A_T_lambda = A_dagger * (Ad * qd + A * (M \ [U_obj; U]));
    f_in = - A_f_t \ A_T_lambda(n_o + 1 : n_o + n_am*3);
    tau_in = A_tau_r \ ( - A_T_lambda(n_o + n_am*3+1 : end) - A_f_r * f_in);

    qdd = M\ ([U_obj ; U] - A_T_lambda);
end

function out = Ry(q)
    out = [cos(q), 0, sin(q);
           0, 1, 0;
           -sin(q), 0, cos(q)];
end