addpath("../dynamics", "../functions", "../../params" )
clear; close all
params = define_params();
m0 = params{1}; I0 = params{2}; mu = params{3}; r= params{4}; d= params{5};
thrust_limit= params{6}; kt = params{7}; c_1 = params{8}; c_2 = params{9}; mass_door = params{10};
handle_factor = params{11}; inertia = params{12}; r_i_ci = params{13}; n = params{14}; dh = params{15}; gravity = params{16};
%%
wn = 1.0; damp = 1.2; % 1 / 1.2
k_p_x = m0 * wn^2;
k_d_x = 2 * damp * sqrt(m0 * k_p_x);  

wn = 3.0; damp = 1.2; % 2 / 1.2
k_p_z = m0 * wn^2;
k_d_z = 2 * damp * sqrt(m0 * k_p_z);  

k_p = diag([k_p_x, k_p_x, k_p_z ]);
k_d = diag([k_d_x, k_d_x, k_d_z ]);

wn = 2; damp = 1.2; % 5 / 1.2
k_R = m0 * wn^2;
k_w = 2 * damp * sqrt(m0 * k_p);  

uncetainty = 1.00; % mass, inertia

dt_sim = 0.001;
N_sim = 15000;

X_hover = [1; 2; 3] * 1e-1; yaw_hover = 0 / 180 *pi; 
[X_des, Xd_des, Xdd_des, yaw_des, yawd_des, yawdd_des] = get_traj_hover(X_hover, yaw_hover, N_sim, dt_sim);
radius = 0.0;  v_z = 0.00;
omega     = 2 * pi * 0.1; 
omega_yaw = 2 * pi * 0.00; % omega_yaw should be small
X_hover = [0.0; 0.0; 0.1]; yaw_hover = 0 / 180 *pi; 
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
B_tau = B_theta(1:3, :);
B_tau_dagger = B_tau' * inv(B_tau * B_tau'); 
B_force = B_theta(4:6, :);
lambda_f_unit = null(B_tau) / norm(null(B_tau));
e_f = B_force * null(B_tau) / norm(B_force * null(B_tau));

e_1 = [1; 0; 0];
e_2 = [0; 1; 0];
e_3 = [0; 0; 1];
%% Passive decomposition
delta_bot = [e_1, e_2];
delta_top = null(delta_bot' * I0);
omega_bot = null(delta_top')';
omega_bot = delta_bot' * I0;
%regularize
delta_bot = inv(I0) * omega_bot' * inv(omega_bot * inv(I0) * omega_bot');
omega_top = inv(delta_top' * I0 * delta_top) * delta_top' * I0;
%omega_top = e_3';

delta = [delta_top delta_bot];
omega = [omega_top; omega_bot];

% Define numerical J
J_num = I0 / norm(I0);

% Define equations
eqs = @(w) [w(2)*(J_num(3,1)*w(1) + J_num(3,2)*w(2) + J_num(3,3)) - ...
             (J_num(2,1)*w(1) + J_num(2,2)*w(2) + J_num(2,3));

             (J_num(1,1)*w(1) + J_num(1,2)*w(2) + J_num(1,3)) - ...
             w(1)*(J_num(3,1)*w(1) + J_num(3,2)*w(2) + J_num(3,3))];

% Solve numerically
w0 = [0; 0];  % Initial guess
options = optimoptions('fsolve', 'Display', 'off');
sol = fsolve(eqs, w0, options);

% Verify
w_des_hat = [sol; 1];
cross_result = cross(w_des_hat, J_num * w_des_hat);
%disp('Solution:'); disp(w_des_hat);
%disp('Cross product:'); disp(cross_result(1:2));  % Should be near [0; 0]

w_des_hat = [0;0;1];
%%
X = [0 ; 0; 0]; Xd = [0 ; 0; 0];
w = [0 ; 0; 0]; wd = [0 ; 0; 0];
R = eye(3, 3);
Rd = R * S(w);
V = [w; R' * Xd];

X_hist = []; Xd_hist = []; w_hist = []; wd_hist = []; R_hist = cell(N_sim, 1); 
w_des_hist = [];wd_des_hist = []; thrusts_hist = [];
v_L_hist = []; v_E_hist = []; v_E_des_hist = [];
e_p_hist = []; e_d_hist = []; e_R_hist = []; e_w_hist = [];
times = [];

for i = 1:N_sim
    %fprintf("time step: %d\n", i);
    % Compute G, C, g
    w = V(1:3);
    v = V(4:6);
    Xd = R * V(4:6);
    G = [I0, zeros(3, 3); zeros(3, 3), m0 * eye(3, 3)];
    C = [S(w) S(v); zeros(3,3) S(w)];
    g = [ zeros(3,1); R' * -m0 * gravity];

    % Geometric controller
    e_p = X - X_des(:, i);
    e_d = Xd - Xd_des(:, i);
    
    w_des = [0 ; 0; yawd_des(:, i)];
    wd_des = [0 ; 0; yawdd_des(:, i)];
    
    w_des = yawd_des(:, i) * w_des_hat;
    wd_des = yawdd_des(:, i) * w_des_hat;
    

    % for the case of e_f = e_3
    b1_des = [cos(yaw_des(i)); sin(yaw_des(i)) ; 0];
    f_des = - m0 * gravity + m0 * Xdd_des(:, i) - k_p * e_p - k_d * e_d;
    b3_des = f_des / norm(f_des);
    b2_des = S(b3_des) * b1_des / norm(S(b3_des) * b1_des);
    R_des = [S(b2_des)*b3_des, b2_des, b3_des];

    e_R = 0.5 * vee(R_des' * R - R' * R_des);
    e_w = w - R' * R_des * w_des;
    tau_des = -k_R * e_R - k_w * e_w + S(w) * I0 * w ...
              -I0 * (S(w) * R' * R_des * w_des - R' * R_des * wd_des);

    lambda_tau = B_tau_dagger * tau_des;
    
    u_4 = (f_des' * R * e_f) / ((B_force * lambda_f_unit)' * e_f);
    lambda_force = u_4 * lambda_f_unit;

    lambda = lambda_tau + lambda_force;

    % Dynamics
    G = G * uncetainty;
    g = g * uncetainty;
    C = C * uncetainty;
    Vd = inv(G) * (B_theta * lambda - g - C * V);
    Rd = R * S(w);
    
    wd = Vd(1:3);
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
    w_des_hist = [w_des_hist, w_des];
    wd_hist = [wd_hist, wd];
    wd_des_hist = [wd_des_hist, wd_des]; 
    R_hist{i} = R;
    thrusts_hist = [thrusts_hist, lambda];
    e_p_hist = [e_p_hist, e_p]; 
    e_d_hist = [e_d_hist, e_d]; 
    e_R_hist = [e_R_hist, e_R]; 
    e_w_hist = [e_w_hist, e_w];

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
for j = 1:4
    plot(times, thrusts_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$f_1$', '$f_2$', '$f_3$', '$f_4$'}, ...
        'Interpreter','latex','FontSize', 12);
title('Thrusts', 'Interpreter', 'latex','FontSize', 14)
grid on
%% Error 
figure('Position',[950 50 800 800]);
colors = lines(6);

%e_p
subplot(2,2,1)
hold on
for j = 1:3
    plot(times, e_p_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_p$', 'Interpreter', 'latex','FontSize', 14)
grid on

%e_d
subplot(2,2,2)
hold on
for j = 1:3
    plot(times, e_d_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_v$', 'Interpreter', 'latex','FontSize', 14)
grid on

%e_R
subplot(2,2,3)
hold on
for j = 1:3
    plot(times, e_R_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_R$', 'Interpreter', 'latex','FontSize', 14)
grid on

%e_w
subplot(2,2,4)
hold on
for j = 1:3
    plot(times, e_w_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_w$', 'Interpreter', 'latex','FontSize', 14)
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