addpath("../dynamics", "../functions", "../../params" )
clear; close all
params = define_params_ver2();
m0 = params{1} / 1.5; I0 = params{2}; mu = params{3}; r = params{4}; d = params{5};
thrust_limit= params{6}; gravity = params{16};
dt_sim = 0.001;
N_sim = 20000;
%% inertia
shape = [1 1 1; % <- x | y
         0 2 0; %      v
         1 1 0];
shape_mass = [m0,  m0, m0;
              m0,  m0 * 1.5, m0;
              m0,  m0, m0];
shape_idx = shape;
[core_row, core_col] = find(shape == 2);
[AMs_rows, AMs_cols] = find(shape ~= 0);

num_AMs = length(AMs_cols);
AM_mass = 0; % mass of shape
AM_inertia = zeros(3, 3); % inertia w.r.t. its com
AM_com = [0; 0; 0];% 0 to com
r_0j = cell(num_AMs, 1); % 0 to j'th module
r_cj = cell(num_AMs, 1); % com to j'th module
I_cj = cell(num_AMs, 1); % inertia of j'th module w.r.t. com of shpe
mass_ams = zeros(num_AMs, 1);
for j =1:length(AMs_cols)
    AM_mass = AM_mass + shape_mass(AMs_rows(j), AMs_cols(j));
    mass_ams(j) = shape_mass(AMs_rows(j), AMs_cols(j));
end
for j =1:length(AMs_cols)
    r_0j{j} = -[(core_col - AMs_cols(j)) *-d ; (core_row - AMs_rows(j)) *d ;0];% p_core to j
    AM_com = AM_com + r_0j{j} * m0/AM_mass;
    shape_idx(AMs_rows(j), AMs_cols(j)) = j;
end
for j =1:length(AMs_cols)
    r_cj{j} = -[(core_col - AMs_cols(j)) *-d ; (core_row - AMs_rows(j)) *d ;0] - AM_com;%p_com to j
    I_cj{j} = I0 + m0 * (r_cj{j}' * r_cj{j} * eye(3,3) - r_cj{j} * r_cj{j}');
    AM_inertia = AM_inertia + I_cj{j};
end
%%
wn = 1; damp = 1.2; % 1, 1.2
kp_M = wn^2; 
kv_M = 2 * damp *sqrt(kp_M);

wn = 1; damp = 1.2; % 1, 1.2
kp_z = wn^2; 
kv_z = 2 * damp *sqrt(kp_z);

kp_M = diag([kp_M, kp_M, kp_z]);
kv_M = diag([kv_M, kv_M, kv_z]);

kw_I = 10; %10

%regulizer
k_lambda = 0.0 / norm(gravity);
k_delta = 1 / norm(gravity);

epsilon = kv_M(1, 1) * 0.3;
alpha = 30 / m0; % 30 / m0
gamma = 50 / m0; % 50 / m0 
beta = diag([1 1 5]) * m0 * norm(gravity) * 3; % 1 1 5 * m0 * norm(gravity) * 3
N_sim_tmp = 10000;

mass_uncertainty = 1.10; 
inertia_uncertainty = 1.10;
thrust_limit = thrust_limit * 3;

%disturbance
sigma = 0.000; mean = 5.0; max_val = 30.0;
disturb_sim = zeros(N_sim, 6);
disturb = mean * [0; 0; 0; 0.5; -0.7; -1.0];
% X, w_estimation error
X_error = zeros(3, num_AMs);
w_error = zeros(3, num_AMs);
sigma_X = 0 / 100; max_X = 0.1; % delicate
sigma_w = 2 / 100; max_w = 0.1; 

rng('shuffle')

X_hover = [1; 2; 3] * 1e-1; yaw_hover = 0 / 180 *pi; 
[X_des, Xd_des, Xdd_des, yaw_des, yawd_des, yawdd_des] = get_traj_hover(X_hover, yaw_hover, N_sim, dt_sim);
radius = 0.3;  v_z = 0.05;
omega     = 2 * pi * 0.1; 
omega_yaw = 2 * pi * 0.05; 
X_hover = [0; 0; 0.5]; yaw_hover = 0 / 180 *pi; 
%[X_des, Xd_des, Xdd_des, yaw_des, yawd_des, yawdd_des] = get_traj_helix(radius, omega, omega_yaw, v_z, X_hover, yaw_hover, N_sim, dt_sim);
%%
B = [r r -r -r;  -r r r -r; mu -mu mu -mu; 1 1 1 1];

e_1 = [1; 0; 0];
e_2 = [0; 1; 0];
e_3 = [0; 0; 1];
%%
X = [0; 0; 0]; Xd = [0; 0; 0]; Xdd = [0; 0; 0];
w = [0; 0; 0]; wd = [0; 0; 0];
R = eye(3, 3);
Rd = R * S(w);
V = [w; R' * Xd];

force_prev = ones(num_AMs) * m0 * 0.7;
delta_hat = zeros(3, num_AMs);
w_des_prev = zeros(3, num_AMs);

X_hist = []; Xd_hist = []; w_hist = []; wd_hist = []; R_hist = cell(N_sim, 1); 
w_des_hist = []; wd_des_hist = []; thrusts_hist = []; F_hat_hist = [];
e_p_hist = []; e_d_hist = [];e_w_hist = []; nu_e_hist = []; delta_hat_hist = []; delta_tilde_hist = [];
times = [];

for i = 1:N_sim_tmp
    % Compute G, C, g
    w = V(1:3);
    v = V(4:6);
    Xd = R * V(4:6);
    G = [AM_inertia, zeros(3, 3); zeros(3, 3), AM_mass * eye(3, 3)];
    C = [S(w) S(v); zeros(3,3) S(w)];
    g = [ zeros(3,1); R' * AM_mass * -gravity];
    
    % Semi-decentralized control
    tau_tot = zeros(3, 1);
    force_tot = 0;
    thrusts = [];
    %estimation error
    X_error_dot = randn(3, num_AMs) * sigma_X;
    X_error = X_error + X_error_dot * dt_sim;
    X_error = min(max(X_error, - max_X), max_X);
    w_error_dot = randn(3, num_AMs) * sigma_w;
    w_error = w_error + w_error_dot * dt_sim;
    w_error = min(max(w_error, - max_w), max_w);

    for j = 1:num_AMs
        
        % position control - backstepping
        mj = mass_ams(j);
        kp_j = mj * kp_M;
        kv_j = mj * kv_M;
        
        e_p  = (1 - X_error(:, j)).*X - X_des(:, i); % could be computed using X_com = X_j - R' * r_cj
        e_pd = (1 - X_error(:, j)).*Xd - Xd_des(:, i);
        e_pdd_hat = gravity + force_prev(j) / mj * R * e_3 - Xdd_des(:, i);
        %e_pdd_hat = Xdd - Xdd_des(:, i);
        
        if i < N_sim
            Xddd_des = (Xdd_des(:, i+1) - Xdd_des(:, i)) / dt_sim;
        end
        
        %regulizer : not stable
        if i > 1
            force_prev(j) = force_prev(j) + k_lambda * (mj / AM_mass * force_tot - force_prev(j));
        end

        nu_ej = force_prev(j) * R *e_3 /mj - Xdd_des(:, i) + kv_j * e_pd / mj ...
                + kp_j * e_p / mj + gravity + delta_hat(:, j) / mj;
        eta = - alpha * mj * nu_ej ...
              - mj * gamma * (e_pd + epsilon * e_p);
        vec = R' * (eta - mj * Xddd_des + kv_j * e_pdd_hat + kp_j * e_pd + delta_hat(:, j));
        w_xj_des = - vec(2) / force_prev(j);
        w_yj_des = vec(1) / force_prev(j);
        forced = vec(3);

        force_j = force_prev(j) + forced * dt_sim;
        force_prev(j) = force_j;
        
        % adaptive control
        delta_hatd = mj / AM_mass * beta * (e_pd + epsilon * e_p + kv_j / mj / gamma * nu_ej); 
        delta_hat(:, j) = delta_hat(:, j) + delta_hatd * dt_sim;

        % rotaion control
        Ij = I_cj{j};
        kw_j = kw_I * norm(Ij);
        
        
        w_des = [w_xj_des; w_yj_des; yawd_des(i)]; 
        if i > 1
            wd_des = w_des - w_des_prev(:, j);
        else
            wd_des = [0; 0; 0];
        end
        w_des_prev(:, j) = w_des;

        wj = (1 - w_error(:, j)).*w;
        e_w = wj - w_des;
        tau_j = S(wj) * Ij * wj + Ij * wd_des - kw_j * e_w;
        
        lambda_j = inv(B) * [tau_j; force_j];
        lambda_j = min(thrust_limit, max(-thrust_limit, lambda_j));
        thrusts = [thrusts; lambda_j];

        tau_j = B(1:3,:) * lambda_j;
        force_j = B(4, :) * lambda_j;
        
        
        tau_tot = tau_tot + tau_j + S(r_cj{j}) * [0; 0; force_j];
        force_tot = force_tot + force_j;
    end

    % generate disturbance
    disturb_dot = randn(6, 1) * sigma;
    disturb = disturb + disturb_dot * dt_sim;
    disturb = min(max(disturb, - max_val), max_val);
    disturb_sim(i,:) = disturb;

    % Dynamics
    G = G * diag([ones(3,1) * inertia_uncertainty; ones(3,1) * mass_uncertainty]);
    g = g * mass_uncertainty;
    Vd = inv(G) * ([tau_tot; 0; 0; force_tot] + disturb - g - C * V);
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
    thrusts_hist = [thrusts_hist, thrusts];

    e_p = X - X_des(:, i); 
    e_pd = Xd - Xd_des(:, i);
    e_w = w - w_des;
   
    e_p_hist = [e_p_hist, e_p]; 
    e_d_hist = [e_d_hist, e_pd]; 
    e_w_hist = [e_w_hist, e_w];
    nu_e_hist = [nu_e_hist, nu_ej];
    delta_hat_hist = [delta_hat_hist, delta_hat(:, j)];
    delta_tilde_hist = [delta_tilde_hist, mj / AM_mass * disturb(4:6) - delta_hat(:, j)];

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
%ylim([-1, 1])
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
%ylim([-1, 1])
grid on

% 7. thrusts
subplot(3,2,5)
hold on
plot(times, thrusts_hist, 'LineWidth', 1.0);
ylim([ -thrust_limit, thrust_limit] )
title('Thrusts', 'Interpreter', 'latex','FontSize', 14)
grid on

%% Error 
figure('Position',[900 50 1000 700]);
colors = lines(6);

%e_p
subplot(3,2,1)
hold on
for j = 1:3
    plot(times, e_p_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_p$', 'Interpreter', 'latex','FontSize', 14)
grid on

%e_d
subplot(3,2,2)
hold on
for j = 1:3
    plot(times, e_d_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_v$', 'Interpreter', 'latex','FontSize', 14)
grid on

%e_w
subplot(3,2,3)
hold on
for j = 1:3
    plot(times, e_w_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_w$', 'Interpreter', 'latex','FontSize', 14)
%ylim([-1, 1])
grid on

%nu_e
subplot(3,2,4)
hold on
for j = 1:3
    plot(times, nu_e_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$\nu_{e,i}$', 'Interpreter', 'latex','FontSize', 14)
grid on

%delta_hat
subplot(3,2,5)
hold on
for j = 1:3
    plot(times, delta_hat_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$\hat{\Delta}_{p,i}$', 'Interpreter', 'latex','FontSize', 14)
grid on

%delta_tilde
subplot(3,2,6)
hold on
for j = 1:3
    plot(times, delta_tilde_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$\tilde{\Delta}_{p,i}$', 'Interpreter', 'latex','FontSize', 14)
grid on