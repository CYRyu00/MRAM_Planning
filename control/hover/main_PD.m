addpath("../dynamics", "../functions", "../../params" )
clear; close all
params = define_params();
m0 = params{1}; I0 = params{2}; mu = params{3}; r= params{4}; d= params{5};
thrust_limit= params{6}; kt = params{7}; c_1 = params{8}; c_2 = params{9}; mass_door = params{10};
handle_factor = params{11}; inertia = params{12}; r_i_ci = params{13}; n = params{14}; dh = params{15}; gravity = params{16};
%%
wn = 2; damp = 0.9; 
kp_pos = m0 * wn^2;
kd_pos = 2 * damp * sqrt(m0 * kp_pos);  

wn = 10; damp = 1.0; 
kp_yaw = m0 * wn^2;
kd_yaw = 2 * damp * sqrt(m0 * kp_yaw);  

K_p = diag([kp_pos kp_pos kp_pos kp_yaw]);
K_d = diag([kd_pos kd_pos kd_pos kd_yaw]);

X_hover = [1; 2; 3] * 1e-2; yaw_hover = 1 / 180 *pi; 
dt_sim = 0.001;
N_sim = 300;
%%
theta = 15 / 180 * pi;
s = sin(theta); c = cos(theta);
B_theta = [mu*s/sqrt(2) - r*c, mu*s/sqrt(2) - r*c, - mu*s/sqrt(2) + r*c, - mu*s/sqrt(2) + r*c;
         - mu*s/sqrt(2) + r*c, mu*s/sqrt(2) - r*c, mu*s/sqrt(2) - r*c, - mu*s/sqrt(2) + r*c;
         mu*c + sqrt(2)*r*s, - mu*c - sqrt(2)*r*s, mu*c + sqrt(2)*r*s, - mu*c - sqrt(2)*r*s;
         s/sqrt(2), - s/sqrt(2), - s/sqrt(2), s/sqrt(2);
         - s/sqrt(2), - s/sqrt(2), s/sqrt(2), s/sqrt(2);
         c, c, c, c];
B_tau = B_theta(1:3, :);
B_force = B_theta(4:6, :);

e_1 = [1; 0; 0];
e_2 = [0; 1; 0];
e_3 = [0; 0; 1];
%%
X = [0 ; 0; 0]; Xd = [0 ; 0; 0];
w = [0 ; 0; 0]; wd = [0 ; 0; 0];
R = eye(3, 3);
Rd = R * S(w);
V = [w; R' * Xd];

X_hist = []; Xd_hist = []; w_hist = []; wd_hist = []; R_hist = cell(N_sim, 1); wd_des_hist = [];
thrusts_hist = [];
v_L_hist = []; v_E_hist = []; v_E_des_hist = [];
times = [];

[X_des, Xd_des, Xdd_des, yaw_des, yawd_des, yawdd_des] = get_traj_hover(X_hover, yaw_hover, N_sim, dt_sim);
%[X_des, Xd_des, Xdd_des, yaw_des, yawd_des, yawdd_des] = get_traj_helix(X_hover, yaw_hover, N_sim, dt_sim);

for i = 1:N_sim
    %fprintf("\n\ntime step: %d\n", i);
    % Compute G, C, g
    w = V(1:3);
    v = V(4:6);
    G = [I0, zeros(3, 3); zeros(3, 3), m0 * eye(3, 3)];
    C = [S(w) S(v); zeros(3,3) S(w)];
    g = [ zeros(3,1); R' * m0 * -gravity];

    % Passive decomposition
    A = [zeros(3,3), R; e_3' * R, zeros(1,3)];
    omega_bot = A;
    delta_top = [R' * e_1, R' * e_2; zeros(3, 2)];
    delta_bot = inv(G) * omega_bot' * inv(omega_bot * inv(G) * omega_bot');
    omega_top = inv(delta_top' * G * delta_top) * delta_top' * G;
    delta = [delta_top, delta_bot];
    
    B_L = delta_top' * B_theta;
    B_E = delta_bot' * B_theta;
    
    v_L = omega_top * V;
    v_E = omega_bot * V;
    
    G_L = delta_top' * G * delta_top;
    G_E = delta_bot' * G * delta_bot;

    Ad =  [zeros(3,3), Rd; e_3' * Rd, zeros(1,3)];
    delta_topd = [Rd' * e_1, Rd' * e_2; zeros(3,2)];
    delta_botd = inv(G) * Ad' * inv(A * inv(G) * A') ...
                 - inv(G) * A' * inv(A * inv(G) * A') ...
                 * (Ad * inv(G) * A' + A * inv(G) * Ad') * inv(A * inv(G) * A');
    deltad = [delta_topd, delta_botd];
    C_pd = delta' * (G * deltad + C * delta);
    C_L  = C_pd(1:2, 1:2);
    C_LE = C_pd(1:2, 3:6);
    C_EL = C_pd(3:6, 1:2);
    C_E  = C_pd(3:6, 3:6);

    g_L = delta_top' * g;
    g_E = delta_bot' * g;
    
    % Desired traj
    yaw = atan2(R(2, 1), R(1, 1));
    h = [X; yaw];
    h_des = [X_des(:, i); yaw_des(:, i)];
    w_des = [0 ; 0; yawd_des(:, i)];
    wd_des = [0 ; 0; yawdd_des(:, i)];
    
    R_des = rpy2rot(0, 0, yaw_des(:, i));
    Rd_des = R_des * S(w_des);

    v_E_des = [Xd_des(:, i); e_3' * R_des * w_des];
    vd_E_des = [Xdd_des(:, i); e_3' * (Rd_des * w_des + R_des * wd_des)];
    
    B_E_lambda = C_EL * v_L + C_E * v_E + g_E ...
                 + G_E * vd_E_des - K_p * (h - h_des) - K_d * (v_E - v_E_des);

    lambda = inv(B_E) * B_E_lambda;

    % Dynamics
    Vd = inv(G) * (B_theta * lambda - g - C * V);
    Rd = R * S(w);
    
    wd = Vd(1:3);
    Xd = R * V(4:6);
    Xdd = R * Vd(4:6);
    
    X  = X + Xd * dt_sim + 0.5 * Xdd * dt_sim^2 ;
    V = V + Vd * dt_sim;
    R  = R + Rd * dt_sim;
    [U_, ~, V_] = svd(R);
    R  = U_ * V_';   
    
    % Log
    X_hist  = [X_hist, X];
    Xd_hist = [Xd_hist, Xd];
    w_hist  = [w_hist, w];
    wd_hist = [wd_hist, wd];
    wd_des_hist = [wd_des_hist, wd_des]; 
    R_hist{i} = R;
    thrusts_hist = [thrusts_hist, lambda];
    v_L_hist = [v_L_hist, v_L];
    v_E_hist = [v_E_hist, v_E];
    v_E_des_hist = [v_E_des_hist, v_E_des];

    times = [times; i * dt_sim];
end
%% State plot
figure('Position',[100 50 1000 900]);
colors = lines(6);

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

% v_L, v_E
subplot(3,2,6)
hold on
for j = 1:2
    plot(times, v_L_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
for j = 1:4
    plot(times, v_E_hist(j, :), 'Color', colors(j + 2,:), 'LineWidth', 1.0);
end
legend({'$v_{L1}$', '$v_{L2}$', '$v_{E1}$', '$v_{E2}$', '$v_{E3}$', '$v_{E4}$'}, ...
        'Interpreter','latex','FontSize', 12);
title('$v_L, v_E$', 'Interpreter', 'latex','FontSize', 14)
grid on

%%
function [Vd, Rd] = FowardDynamics(R, V, m, J, gravity, B, lambda)
    w = V(1:3);
    v = V(4:6);
    G = [J, zeros(3, 3); zeros(3, 3), m * eye(3, 3)];
    C = [S(w) S(v); zeros(3,3) S(w)];
    g = R' * m * gravity;

    Vd = inv(G) * (B * lambda - g - C * V);
    Rd = R * S(w);
end