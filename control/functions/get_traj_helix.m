function [X_des, Xd_des, Xdd_des, ...
          yaw_des, yawd_des, yawdd_des] = get_traj_helix(X_hover, yaw_hover, ...
                                                         N_sim, dt_sim)
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

    %% Trajectory design parameters (수정 가능)
    r     = 0.5;        % [m]   수평 반지름
    omega = 0.6;        % [rad/s]  회전 각속도 (2π rad ≈ 한 바퀴)
    v_z   = 0.05;       % [m/s]  수직 상승(음수면 하강)

    %% 사전 할당 (pre-allocation)
    X_des   = zeros(3, N_sim);
    Xd_des  = zeros(3, N_sim);
    Xdd_des = zeros(3, N_sim);

    yaw_des   = zeros(1, N_sim);
    yawd_des  = zeros(1, N_sim);
    yawdd_des = zeros(1, N_sim);

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

        %── yaw  (진행 방향을 바라보도록 설정) ───────────────────
        %   yaw = atan2(ŷ속도, ẋ속도)  ;  필요시 yaw_hover 오프셋 추가
        yaw   = atan2(ydot, xdot) + yaw_hover;
        yawd  = omega;          % 등속 회전이므로 상수
        yawdd = 0;

        yaw_des(k)   = wrapToPi(yaw);
        yawd_des(k)  = yawd;
        yawdd_des(k) = yawdd;
    end
end
