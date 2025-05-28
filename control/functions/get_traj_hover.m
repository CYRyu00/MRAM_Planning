function [X_des, Xd_des, Xdd_des, yaw_des, yawd_des, yawdd_des] = get_traj_hover(X_hover, yaw_hover, N_sim, dt_sim)
    X_des = []; Xd_des = []; Xdd_des = []; 
    yaw_des = []; yawd_des = []; yawdd_des = [];

    X_prev = X_hover; Xd_prev = zeros(3,1);
    yaw_prev = yaw_hover; yawd_prev = 0;
    for i = 1:N_sim
        X = X_hover;
        Xd = (X - X_prev) / dt_sim; 
        Xdd = (Xd - Xd_prev) / dt_sim; 
        X_des = [X_des, X];
        Xd_des = [Xd_des, Xd];
        Xdd_des = [Xdd_des, Xdd];
        X_prev = X; Xd_prev = Xd;

        yaw = yaw_hover;
        yawd = (yaw - yaw_prev) / dt_sim;
        yawdd = (yawd - yawd_prev) / dt_sim;
        yaw_des = [yaw_des, yaw];
        yawd_des = [yawd_des, yawd];
        yawdd_des = [yawdd_des, yawdd];
        yaw_prev = yaw;
        yawd_prev = yawd;
    end
end