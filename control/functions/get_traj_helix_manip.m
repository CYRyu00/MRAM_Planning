function  [X_des, Xd_des, Xdd_des, ...
           R_e_des, w_e_des] = get_traj_helix_manip(r, omega, v_z, rpyd, ...
                                                    X_hover, rpy_hover, N_sim, dt_sim)
% GET_TRAJ_HELIX  Helical reference trajectory (position–velocity–acceleration & yaw)
%
%   INPUTS
%     X_hover  – 3×1 [m]   : 중심(시작) 위치 (x, y, z)
%     yaw_hover– 1×1 [rad] : 시작 yaw角 (0 rad = x-축 정방향)
%     N_sim    – 1×1       : 시뮬레이션 스텝 수
%     dt_sim   – 1×1 [s]   : 스텝 시간
%
%   OUTPUTS
%     X_des    – 3×N_sim   : 위치
%     Xd_des   – 3×N_sim   : 속도
%     Xdd_des  – 3×N_sim   : 가속도
%     yaw_des  – 1×N_sim   : yaw
%     yawd_des – 1×N_sim   : yaw 속도
%     yawdd_des– 1×N_sim   : yaw 가속도

    %% 사전 할당 (pre-allocation)
    X_des   = zeros(3, N_sim);
    Xd_des  = zeros(3, N_sim);
    Xdd_des = zeros(3, N_sim);

    R_e_des = cell(N_sim, 1);
    w_e_des = [];
    R_e_prev = rpy2rot(rpy_hover(1), rpy_hover(2), rpy_hover(3));

    %% 메인 루프
    for k = 1:N_sim
        t = (k-1) * dt_sim;     % 현재 시간 [s]
        theta = omega * t;      % 진행 각도 [rad]

        %── 위치 ────────────────────────────────────────────────
        x = X_hover(1) + r * cos(theta);
        y = X_hover(2) + r * sin(theta);
        z = X_hover(3) + v_z * t;
        X_des(:,k) = [x; y; z];

        %── 속도 (1차 미분) ──────────────────────────────────────
        xdot = -r * omega * sin(theta);
        ydot =  r * omega * cos(theta);
        zdot =  v_z;
        Xd_des(:,k) = [xdot; ydot; zdot];

        %── 가속도 (2차 미분) ────────────────────────────────────
        xddot = -r * omega^2 * cos(theta);
        yddot = -r * omega^2 * sin(theta);
        zddot =  0;
        Xdd_des(:,k) = [xddot; yddot; zddot];

        roll = rpy_hover(1) + dt_sim * (k-1) * rpyd(1);
        pitch = rpy_hover(2) + dt_sim * (k-1) * rpyd(1);
        yaw = rpy_hover(3) + dt_sim * (k-1) * rpyd(1);
        R_e = rpy2rot(roll, pitch, yaw);
        Rd = (R_e - R_e_prev) / dt_sim;
        w_e = R_e' * Rd;
        R_e_prev = R_e;
        R_e_des{k} = R_e;
        w_e_des = [w_e_des, w_e];
    end
end
