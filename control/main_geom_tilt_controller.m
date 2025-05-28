addpath("dynamics", "functions", "../params" )
clear; close all
params = define_params();
m0 = params{1}; I0 = params{2}; mu = params{3}; r= params{4}; d= params{5};
thrust_limit= params{6}; kt = params{7}; c_1 = params{8}; c_2 = params{9}; mass_door = params{10};
handle_factor = params{11}; inertia = params{12}; r_i_ci = params{13}; n = params{14}; dh = params{15}; gravity = params{16};
%%
%m0 = 4.34; I0 = diag([0.082, 0.084, 0.130]);
wn = 2; damp = 0.8; % 30: 2/ 0.7
kp = m0 * wn^2;
kv = 2 * damp * sqrt(m0*kp);  

wn = 7; damp = 0.5;% 30:7 / 0.5
kR = max(max(I0)) * wn^2;
kw = 2 * damp * sqrt(max(max(I0)) * kR);  

%kp = 9; kv = 15;
%kR = 5; kw = 0.8;

fprintf("\nkp: %f, kv %f\n", kp, kv);
fprintf("\nkR: %f, kw %f\n", kR, kw);
X_hover = [1; 2; 3] * 1e-1; yaw_hover = 10 / 180 *pi; 
L_f = eye(3) * 10;
L_tau = eye(3) * 100;

theta = 30 / 180 * pi;
thrust_max = thrust_limit * 1e2;
thrust_min = thrust_limit * -1e2;
%%
N = 100; dt = 0.1;
dt_sim = 0.01;
N_sim = N * dt / dt_sim;
times = [];

s = sin(theta); c = cos(theta);
A_theta = [mu*s/sqrt(2) - r*c, mu*s/sqrt(2) - r*c, - mu*s/sqrt(2) + r*c, - mu*s/sqrt(2) + r*c;
         - mu*s/sqrt(2) + r*c, mu*s/sqrt(2) - r*c, mu*s/sqrt(2) - r*c, - mu*s/sqrt(2) + r*c;
         mu*c + sqrt(2)*r*s, - mu*c - sqrt(2)*r*s, mu*c + sqrt(2)*r*s, - mu*c - sqrt(2)*r*s;
         s/sqrt(2), - s/sqrt(2), - s/sqrt(2), s/sqrt(2);
         - s/sqrt(2), - s/sqrt(2), s/sqrt(2), s/sqrt(2);
         c, c, c, c];
A_tau = A_theta(1:3, :);
A_force = A_theta(4:6, :);

X = [0 ; 0; 0]; Xd = [0 ; 0; 0];
w = [0 ; 0; 0]; wd = [0 ; 0; 0];
R = eye(3, 3);
thrusts = ones(4, 1) * m0 * norm(gravity) /4;
f_hat = zeros(3, 1); tau_hat = zeros(3, 1); int_wd_des = zeros(3, 1);

X_hist = []; Xd_hist = []; w_hist = []; wd_hist = []; R_hist = cell(N_sim, 1); wd_des_hist = [];
f_hat_hist = []; tau_hat_hist = []; f_e_hist = []; tau_e_hist = [];
thrusts_hist = [];
ep_hist = []; ev_hist = []; eR_hist = []; ew_hist = [];

%[X_des, Xd_des, Xdd_des, yaw_des, yawd_des, yawdd_des] = get_traj_hover(X_hover, yaw_hover, N_sim, dt_sim);
[X_des, Xd_des, Xdd_des, yaw_des, yawd_des, yawdd_des] = get_traj_helix(X_hover, yaw_hover, N_sim, dt_sim);

for i = 1:N_sim
    fprintf("\n\ntime step: %d\n", i);
    % External wrench & wrench estimator
    % TODO estimation by observed w_dot, w or noisy
    f_e = 0 * [0; 0; 0.02] * sin(i/100); 
    tau_e = 0 * [3; 2; 1] * 1e-4 * sin(i/100);
    f_hat  = f_hat  + L_f * (f_e  - f_hat) * dt_sim;
    tau_hat = tau_hat + L_tau * (tau_e - tau_hat) * dt_sim;
    
    f_hat_hist   = [f_hat_hist, f_hat];
    tau_hat_hist = [tau_hat_hist, tau_hat];
    f_e_hist     = [f_e_hist, f_e];
    tau_e_hist   = [tau_e_hist, tau_e];
    
    % Geometric Controller w/ tilted quadrotor
    % Position control
    ep = X - X_des(:, i);
    ev = Xd - Xd_des(:, i);
    f_des = - kp * ep - kv * ev + m0 * norm(gravity) * [0; 0; 1] + m0 * Xdd_des(:, i) - f_hat;

    eta = (A_force * thrusts) / norm(A_force * thrusts);
    force = f_des' * R * eta; 
    
    fprintf("force: %f\n", force );
    fprintf("thrust_prev : %f, %f, %f, %f\n", thrusts(1), thrusts(2), thrusts(3), thrusts(4))
    fprintf("eta : %f, %f, %f\n", eta(1), eta(2), eta(3))
    % Compute R_des
    rho = asin( norm(cross(eta, f_des)) / norm(f_des));
    zeta = cross(eta, f_des) / norm(cross(eta, f_des));
    
    % TODO: how to get beta    
    R_des = get_R_des(yaw_des(i), rho, zeta, eta); 

    % Attitude Control
    eR = 0.5 * vee(R_des' * R - R' * R_des);
    w_des = yawd_des(:, i) * eta;
    wd_des = yawdd_des(:, i) * eta;
    ew = w - R' * R_des * w_des;
    tau_des = - kR * eR - kw * ew + S(w) * I0 * w ...
              - I0 * (S(w) * R' * R_des * w_des - R' * R_des * wd_des) - tau_hat;
    
    % Dynamics
    fprintf("tau_des: %f,  %f,  %f\n", tau_des(1), tau_des(2), tau_des(3));
    thrusts = solve_QP_mapping(A_force, A_tau, tau_des, force, thrust_min, thrust_max, R, eta);
    fprintf("thrust : %f, %f, %f, %f\n", thrusts(1), thrusts(2), thrusts(3), thrusts(4))
    thrusts_hist = [thrusts_hist, thrusts];

    [Xdd, wd_fd, Rd] = ForwardDynamics(X, Xd, w, wd, R, f_e, tau_e, thrusts, A_theta, params);

    X  = X + Xd * dt_sim + 0.5 * Xdd * dt_sim^2 ;
    Xd = Xd + Xdd * dt_sim;
    w  = w + wd_fd * dt_sim;
    wd = wd_fd;
    R  = R + Rd * dt_sim;
    [U, ~, V] = svd(R);
    R  = U * V';   
    
    X_hist  = [X_hist, X];
    Xd_hist = [Xd_hist, Xd];
    w_hist  = [w_hist, wd];
    wd_hist = [wd_hist, wd];
    wd_des_hist = [wd_des_hist, wd_des]; 
    R_hist{i} = R;
    ep_hist = [ep_hist, ep];
    ev_hist = [ev_hist, ev];
    eR_hist = [eR_hist, eR];
    ew_hist = [ew_hist, ew];
    times = [times; i * dt_sim];
end

%% State plot
figure('Position',[100 50 800 900]);
colors = lines(4);

subplot(3,2,1)
hold on
for j = 1:3
    plot(times, X_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    plot(times, X_des(j, 1:i), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
end
legend({'$X_x$', '$X_x^{\mathrm{des}}$', ...
        '$X_y$', '$X_y^{\mathrm{des}}$', ...
        '$X_z$', '$X_z^{\mathrm{des}}$'}, ...
        'Interpreter','latex','FontSize', 12);
title('$\mathbf{X}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
grid on
% 2. Xd plot
subplot(3,2,2)
hold on
for j = 1:3
    plot(times, Xd_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    plot(times, Xd_des(j, 1:i), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
end
legend({'$\dot{X}_x$', '$\dot{X}_x^{\mathrm{des}}$', ...
        '$\dot{X}_y$', '$\dot{X}_y^{\mathrm{des}}$', ...
        '$\dot{X}_z$', '$\dot{X}_z^{\mathrm{des}}$'}, ...
        'Interpreter','latex','FontSize', 12);
title('$\dot{\mathbf{X}}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
grid on
% 3. w plot
subplot(3,2,3)
hold on
for j = 1:3
    plot(times, w_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
plot(times, yawd_des(1, 1:i), '--', 'Color', colors(3,:), 'LineWidth', 2.5);  % w_z des
legend({'$\omega_x$', '$\omega_y$', '$\omega_z$', '$\omega_z^{\mathrm{des}}$'}, ...
       'Interpreter', 'latex','FontSize', 12);
title('${\omega}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
grid on
% 4. wd plot
subplot(3,2,4)
hold on
for j = 1:3
    plot(times, wd_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    plot(times, wd_des_hist(j, 1:i), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
end
legend({'$\dot{\omega}_x$', '$\dot{\omega}_x^{\mathrm{des}}$', ...
        '$\dot{\omega}_y$', '$\dot{\omega}_y^{\mathrm{des}}$', ...
        '$\dot{\omega}_z$', '$\dot{\omega}_z^{\mathrm{des}}$'}, ...
        'Interpreter','latex','FontSize', 12);
title('$\dot{{\omega}}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
grid on

% 7. thrusts
subplot(3,2,5)
hold on
for j = 1:4
    plot(times, thrusts_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$f_1$', '$f_2$', '$f_3$', '$f_4$'}, ...
        'Interpreter','latex','FontSize', 12);
title('Thrusts', 'Interpreter', 'latex','FontSize', 14)
grid on
%% Error plot
figure('Position',[900 100 800 800]);
colors = lines(3);

subplot(3,2,1);  hold on;
for j = 1:3
    plot(times, ep_hist(j,:), 'Color', colors(j,:), 'LineWidth', 1.3);
end
legend({'$e_{p_x}$','$e_{p_y}$','$e_{p_z}$'}, ...
       'Interpreter','latex','FontSize', 12);
ylabel('[m]','Interpreter','latex');
title('Position Error $\mathbf e_p$', 'Interpreter','latex','FontSize', 14);
grid on

subplot(3,2,2);  hold on;
for j = 1:3
    plot(times, ev_hist(j,:), 'Color', colors(j,:), 'LineWidth', 1.3);
end
legend({'$e_{v_x}$','$e_{v_y}$','$e_{v_z}$'}, ...
       'Interpreter','latex','FontSize', 12);
ylabel('[$\mathrm{m/s}$]','Interpreter','latex');
title('Velocity Error $\mathbf e_{v}$', 'Interpreter','latex','FontSize', 14);
grid on

subplot(3,2,3);  hold on;
for j = 1:3
    plot(times, eR_hist(j,:), 'Color', colors(j,:), 'LineWidth', 1.3);
end
legend({'$e_{R_x}$','$e_{R_y}$','$e_{R_z}$'}, ...
       'Interpreter','latex','FontSize', 12);
ylabel('[$\mathrm{m/s}$]','Interpreter','latex');
title('Attitude Error $\mathbf e_{R}$', 'Interpreter','latex','FontSize', 14);
grid on

subplot(3,2,4);  hold on;
for j = 1:3
    plot(times, ew_hist(j,:), 'Color', colors(j,:), 'LineWidth', 1.3);
end
legend({'$e_{\omega_x}$','$e_{\omega_y}$','$e_{\omega_z}$'}, ...
       'Interpreter','latex','FontSize', 12);
ylabel('[$\mathrm{m/s}$]','Interpreter','latex');
title('Angular Velocity Error $\mathbf e_{\omega}$', 'Interpreter','latex','FontSize', 14);
grid on

% force estimation error
f_err = f_e_hist - f_hat_hist;
tau_err = tau_e_hist - tau_hat_hist;

subplot(3,2,5);  hold on;
for j = 1:3
    plot(times, f_err(j,:), 'Color', colors(j,:), 'LineWidth', 1.3);
end
legend({'$e_{f_x}$','$e_{f_y}$','$e_{f_z}$'}, ...
       'Interpreter','latex','FontSize', 12);
ylabel('[$\mathrm{N}$]','Interpreter','latex');
title('Force estimation Error $\mathbf e_{f}$', ...
      'Interpreter','latex','FontSize', 14);
grid on
% moment estimation error
subplot(3,2,6);  hold on;
for j = 1:3
    plot(times, tau_err(j,:), 'Color', colors(j,:), 'LineWidth', 1.3);
end
legend({'$e_{\tau_x}$','$e_{\tau_y}$','$e_{\tau_z}$'}, ...
       'Interpreter','latex','FontSize', 14);
ylabel('[$\mathrm{N}$]','Interpreter','latex');
title('Moment estimation Error $\mathbf e_{f}$', ...
      'Interpreter','latex','FontSize', 14);
grid on

sgtitle('Tracking Errors','Interpreter','latex','FontSize', 16);
%%
function [Xdd, wd_fd, Rd] = ForwardDynamics(X, Xd, w, wd, R, f_e, tau_e, thrusts, A_theta, params)
    m0 = params{1}; I0 = params{2}; mu = params{3}; r= params{4}; gravity = params{16};
    %m0 = 4.34; I0 = diag([0.082, 0.084, 0.130]);

    wrench = A_theta * thrusts;
    tau = wrench(1:3);
    force = wrench(4:6);
    fprintf("Real tau; force : %f, %f, %f, %f\n", tau(1), tau(2), tau(3), norm(force));
    Xdd = (R * force - m0 * norm(gravity) * [0; 0; 1] + f_e) / m0;
    wd_fd = inv(I0) * (tau + tau_e - S(w) * I0 * w);
    Rd = R * S(w);
end

function R_des = get_R_des(yaw_des, rho, zeta, eta)
    R1 = Rot_ang_axis(rho, zeta);
    b_psi = [cos(yaw_des); sin(yaw_des); 0];
    E2 = diag([1 1 0]);
    v1 = E2 * b_psi / norm(E2 * b_psi);
    diff_min = inf;

    beta_arr = (-1:0.01:1) * pi;
    for beta_cand = beta_arr
        R_des_cand = R1 * Rot_ang_axis(beta_cand, eta);
        v2 = E2 * R_des_cand * [1; 0 ; 0]; v2 = v2 / norm(v2);
        diff = norm(v1 - v2);
        if diff < diff_min
            R_des = R_des_cand;
            diff_min = diff;
        end
    end
    fprintf("diff_min of beta : %f\n", diff_min)
end