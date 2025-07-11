addpath("../dynamics", "../functions", "../../params" )
clear; close all
params = define_params();
m0 = params{1}; I0 = params{2}; mu = params{3}; r= params{4}; d= params{5};
thrust_limit= params{6}; kt = params{7}; c_1 = params{8}; c_2 = params{9}; mass_door = params{10};
handle_factor = params{11}; inertia = params{12}; r_i_ci = params{13}; n = params{14}; dh = params{15}; gravity = params{16};
m_duo = 2 * m0;
D = [0.5 * d; 0; 0];
I_duo = 2 * I0 + 2 * m0 * (D' * D * eye(3,3) - D * D');
dt_sim = 0.001;
N_sim = 15000;
%%
wn = 4.0; damp = 1.2; % 4 / 1.2
k_p_x = m_duo * wn^2;
k_d_x = 2 * damp * sqrt(m_duo * k_p_x);  

wn = 2.0; damp = 1.2; % 2 / 1.2
k_p_z = m_duo * wn^2;
k_d_z = 2 * damp * sqrt(m_duo * k_p_z);  

k_p = diag([k_p_x, k_p_x, k_p_z]);
k_d = diag([k_d_x, k_d_x, k_d_z]);

wn = 1.2; damp = 0.9; % 1.2 / 0.9
k_R = m_duo * wn^2;
k_w = 2 * damp * sqrt(m_duo * k_R);  
scale = 5.0; gamma = 0.1; %0.5 / 0.1
k_I = k_w *scale;

K_o_r = [norm(I_duo), norm(I_duo), norm(I_duo)];
K_o_p = [m_duo, m_duo, m_duo] *1e1; %1e1
K_o = diag([K_o_r, K_o_p]) * 1e0; %1e0

uncetainty = 1.10; % mass, inertia
sigma = 0.1; mean = 0.10; max_val = 0.3;
disturb_sim = zeros(N_sim, 6);
disturb = mean * [0; 0; 0; 0; 0; 1];
rng('shuffle')

X_hover = [1; 2; 3] * 1e0; yaw_hover = 100 / 180 *pi; 
[X_des, Xd_des, Xdd_des, yaw_des, yawd_des, yawdd_des] = get_traj_hover(X_hover, yaw_hover, N_sim, dt_sim);
radius = 0.2;  v_z = 0.05;
omega     = 2 * pi * 0.1; 
omega_yaw = 2 * pi * 0.05; 
X_hover = [0; 0; 1]; yaw_hover = 0 / 180 *pi; 
[X_des, Xd_des, Xdd_des, yaw_des, yawd_des, yawdd_des] = get_traj_helix(radius, omega, omega_yaw, v_z, X_hover, yaw_hover, N_sim, dt_sim);
%%
theta = 15 / 180 * pi;
s = sin(theta); c = cos(theta);
B_theta = [mu*s/sqrt(2) - r*c, mu*s/sqrt(2) - r*c, - mu*s/sqrt(2) + r*c, - mu*s/sqrt(2) + r*c;
     - mu*s/sqrt(2) + r*c, mu*s/sqrt(2) - r*c, mu*s/sqrt(2) - r*c, - mu*s/sqrt(2) + r*c;
     mu*c + sqrt(2)*r*s, - mu*c - sqrt(2)*r*s, mu*c + sqrt(2)*r*s, - mu*c - sqrt(2)*r*s;
     s/sqrt(2), - s/sqrt(2), - s/sqrt(2), s/sqrt(2);
     - s/sqrt(2), - s/sqrt(2), s/sqrt(2), s/sqrt(2);
     c, c, c, c];

p1 = [-0.5 * d; 0; 0]; p2 = [0.5 * d; 0; 0]; % p(:,3) = [D; 0; 0];
R1 = eye(3,3); R2 = [0 -1 0;1 0 0; 0 0 1]; % R{3} = eye(3,3);

B_duo = [Ad(R1, p1) * B_theta, Ad(R2, p2) * B_theta];
B_dagger = B_duo' * inv(B_duo * B_duo'); % least square

m_duo = 2 * m0;
r = [0.5 * d; 0; 0];
I_duo = 2 * I0 + 2 * m0 * (r' * r * eye(3,3) - r * r');

e_1 = [1; 0; 0];
e_2 = [0; 1; 0];
e_3 = [0; 0; 1];
%%
X = [0 ; 0; 0]; Xd = [0 ; 0; 0];
w = [0 ; 0; 0]; wd = [0 ; 0; 0];
R = eye(3, 3);
Rd = R * S(w);
V = [w; R' * Xd];

G = [I_duo, zeros(3, 3); zeros(3, 3), m_duo * eye(3, 3)];
p_0 = G * V;
F_prev = zeros(6,1);
F_hat = zeros(6,1);
integral = zeros(6,1);
e_I = zeros(3,1);

X_hist = []; Xd_hist = []; w_hist = []; wd_hist = []; R_hist = cell(N_sim, 1); 
w_des_hist = [];wd_des_hist = []; thrusts_hist = []; F_hat_hist = [];
e_p_hist = []; e_d_hist = []; e_R_hist = []; e_w_hist = []; e_I_hist =[]; e_F_hist = [];
times = [];

for i = 1:N_sim
    % fprintf("time step: %d\n", i);
    % generate disturbance
    disturb_dot = randn(6, 1) * sigma;
    disturb = disturb + disturb_dot * dt_sim;
    disturb = min(max(disturb, - max_val), max_val);
    disturb_sim(i,:) = disturb;

    % Compute G, C, g
    w = V(1:3);
    v = V(4:6);
    Xd = R * V(4:6);
    G = [I_duo, zeros(3, 3); zeros(3, 3), m_duo * eye(3, 3)];
    C = [S(w) S(v); zeros(3,3) S(w)];
    g = [ zeros(3,1); R' * m_duo * -gravity];

    % Impedence controller
    e_p = X - X_des(:, i);
    e_d = Xd - Xd_des(:, i);
    
    force_ = R' * ( -m_duo * gravity + m_duo * Xdd_des(:, i) - k_d * e_d - k_p * e_p);

    w_des = [0 ; 0; yawd_des(:, i)];
    wd_des = [0 ; 0; yawdd_des(:, i)];    

    R_des = rpy2rot(0, 0, yaw_des(i));

    e_R = 0.5 * vee(R_des' * R - R' * R_des);
    e_w = w - R' * R_des * w_des;
    e_I = e_I + (e_w + gamma * e_R) * dt_sim;
    tau_ = -k_R * e_R - k_w * e_w - k_I * e_I + S(w) * I_duo * w ...
           -I_duo * (S(w) * R' * R_des * w_des - R' * R_des * wd_des);
    %MBO
    p = G * V;
    integral = integral + (F_prev + F_hat + C' * V - g) * dt_sim;
    F_hat = K_o * (p - p_0 - integral);
    
    e_F = F_hat - disturb;

    tau_ = tau_ - F_hat(1:3);
    force_ = force_ - F_hat(4:6);

    lambda = B_dagger * [tau_; force_];

    F_prev = [tau_; force_];

    % Dynamics
    G = G * uncetainty;
    g = g * uncetainty;
    C = C * uncetainty;
    Vd = inv(G) * (B_duo * lambda + disturb - g - C * V);
    Rd = R * S(w);
    
    wd = Vd(1:3);
    Xdd = R * Vd(4:6);
    
    X = X + Xd * dt_sim + 0.5 * Xdd * dt_sim^2 ;
    V = V + Vd * dt_sim;
    R = R + Rd * dt_sim;
    [U_, ~, V_] = svd(R);
    R  = U_ * V_';
    
    % Log
    X_hist  = [X_hist, X];
    Xd_hist = [Xd_hist, Xd];
    w_hist  = [w_hist, w];
    w_des_hist = [w_des_hist, w_des];
    wd_hist = [wd_hist, wd];
    wd_des_hist = [wd_des_hist, wd_des]; 
    R_hist{i} = R;
    thrusts_hist = [thrusts_hist, lambda];
    F_hat_hist = [F_hat_hist , F_hat];
    e_p_hist = [e_p_hist, e_p]; 
    e_d_hist = [e_d_hist, e_d]; 
    e_R_hist = [e_R_hist, e_R]; 
    e_w_hist = [e_w_hist, e_w];
    e_I_hist = [e_I_hist, e_I];
    e_F_hist = [e_F_hist, e_F];

    times = [times; i * dt_sim];
end
%% State plot
figure('Position',[100 50 800 900]);
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
    plot(times, w_des_hist(j, 1:i), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
end
legend({'$\omega_x$', '$\omega_x^{\mathrm{des}}$', ...
        '$\omega_y$', '$\omega_y^{\mathrm{des}}$', ...
        '$\omega_z$', '$\omega_z^{\mathrm{des}}$'}, ...
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
plot(times, thrusts_hist, 'LineWidth', 1.0);
title('Thrusts', 'Interpreter', 'latex','FontSize', 14)
grid on

% 8. F_hat
subplot(3,2,6)
hold on
plot(times, F_hat_hist, 'LineWidth', 1.0);
legend({'$\tau_x$', '$\tau_y$', '$\tau_z$', ...
        '$f_x$', '$f_y$', '$f_z$'}, ...
        'Interpreter','latex','FontSize', 12);
title('$\hat{F}$', 'Interpreter', 'latex','FontSize', 14)
grid on
%% Error 
figure('Position',[900 50 1000 700]);
colors = lines(6);

%e_p
subplot(2,3,1)
hold on
for j = 1:3
    plot(times, e_p_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_p$', 'Interpreter', 'latex','FontSize', 14)
grid on

%e_d
subplot(2,3,2)
hold on
for j = 1:3
    plot(times, e_d_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_v$', 'Interpreter', 'latex','FontSize', 14)
grid on

%e_F
subplot(2,3,3)
hold on
plot(times, e_F_hist, 'LineWidth', 1.0);
legend({'$\tau_x$', '$\tau_y$', '$\tau_z$', ...
        '$f_x$', '$f_y$', '$f_z$'}, ...
        'Interpreter','latex','FontSize', 12);
title('$e_F$', 'Interpreter', 'latex','FontSize', 14)
grid on

%e_R
subplot(2,3,4)
hold on
for j = 1:3
    plot(times, e_R_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_R$', 'Interpreter', 'latex','FontSize', 14)
grid on

%e_w
subplot(2,3,5)
hold on
for j = 1:3
    plot(times, e_w_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_w$', 'Interpreter', 'latex','FontSize', 14)
grid on

%e_I
subplot(2,3,6)
hold on
for j = 1:3
    plot(times, e_I_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_I$', 'Interpreter', 'latex','FontSize', 14)
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