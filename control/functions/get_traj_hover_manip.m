function [X_des, Xd_des, Xdd_des, Xddd_des, R_e_des, w_e_des, wd_e_des] = get_traj_hover_manip(X_hover, rpy_hover, N_sim, dt_sim)
    X_des = []; Xd_des = []; Xdd_des = []; Xddd_des = []; 
    R_e_des = cell(N_sim, 1); w_e_des = []; wd_e_des = [];

    X_prev = X_hover; Xd_prev = zeros(3,1); 
    R_e_prev = rpy2rot(rpy_hover(1), rpy_hover(2), rpy_hover(3));
    for i = 1:N_sim
        X = X_hover;
        Xd = (X - X_prev) / dt_sim; 
        Xdd = (Xd - Xd_prev) / dt_sim; 
        X_des = [X_des, X];
        Xd_des = [Xd_des, Xd];
        Xdd_des = [Xdd_des, Xdd];
        Xddd_des = [Xddd_des, [0;0;0]];
        X_prev = X; Xd_prev = Xd;

        R_e = rpy2rot(rpy_hover(1), rpy_hover(2), rpy_hover(3));
        Rd = (R_e - R_e_prev) / dt_sim;
        w_e = vee(R_e' * Rd);
        R_e_prev = R_e;
        R_e_des{i} = R_e;
        w_e_des = [w_e_des, w_e];
        wd_e_des = [wd_e_des, [0;0;0]];
    end
end