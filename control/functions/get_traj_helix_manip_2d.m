function  [X_des, Xd_des, Xdd_des, Xddd_des, ...
           R_e_des, w_e_des, wd_e_des] = get_traj_helix_manip_2d(r, omega, v_z, rpyd, ...
                                                    X_hover, rpy_hover, N_sim, dt_sim)
    X_des   = zeros(3, N_sim);
    Xd_des  = zeros(3, N_sim);
    Xdd_des = zeros(3, N_sim);
    Xddd_des = zeros(3, N_sim);

    R_e_des = cell(N_sim, 1);
    w_e_des = [];
    wd_e_des = [];
    R_e_prev = rpy2rot(rpy_hover(1), rpy_hover(2), rpy_hover(3));
    w_e_prev = [0; 0; 0];
    roll = rpy_hover(1);
    pitch = rpy_hover(2);
    yaw = rpy_hover(3);

    for k = 1:N_sim
        t = (k-1) * dt_sim;     % 현재 시간 [s]
        theta = omega * t;      % 진행 각도 [rad]

        %── 위치 ────────────────────────────────────────────────
        x = X_hover(1) + r * sin(theta);
        y = X_hover(2);
        z = X_hover(3) + v_z * t;
        X_des(:,k) = [x; y; z];

        %── 속도 (1차 미분) ──────────────────────────────────────
        xdot = r * omega * cos(theta);
        ydot = 0;
        zdot = v_z;
        Xd_des(:,k) = [xdot; ydot; zdot];

        %── 가속도 (2차 미분) ────────────────────────────────────
        xddot = -r * omega^2 * sin(theta);
        yddot = 0;
        zddot =  0;
        Xdd_des(:,k) = [xddot; yddot; zddot];
        
        % jerk
        xdddot = -r * omega^3 * cos(theta);
        ydddot = 0;
        zdddot = 0;
        Xddd_des(:,k) = [xdddot; ydddot; zdddot];
        
        % rotation
        R_e = rpy2rot(roll, pitch, yaw);
        Rd = (R_e - R_e_prev) / dt_sim;
        w_e = vee(R_e' * Rd);
        wd_e = (w_e - w_e_prev) / dt_sim;
        
        w_e_prev = w_e;
        R_e_prev = R_e;

        roll = roll + rpyd(1) * dt_sim;
        pitch = pitch + rpyd(2) * dt_sim;
        yaw = yaw + rpyd(3) * dt_sim;

        if roll > pi/4 rpyd(1) = rpyd(1) * (1 - dt_sim/2); end   
        if pitch > pi/4 rpyd(2) = rpyd(2) *(1 - dt_sim/2); end
        if yaw > pi/4 rpyd(3) = rpyd(3) * (1 - dt_sim/2); end

        R_e_des{k} = R_e;
        w_e_des = [w_e_des, w_e];
        wd_e_des = [wd_e_des, wd_e];
    end
end
