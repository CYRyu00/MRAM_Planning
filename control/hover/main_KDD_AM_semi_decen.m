addpath("../dynamics", "../functions", "../../params" )
clear; close all
params = define_params_ver2();
mq = params{1}; Iq = params{2}; mu = params{3}; r = params{4}; d = params{5};
thrust_limit= params{6}; gravity = params{16};
mb = params{17}; cb = params{18}; Ib = params{19};
mt = params{20}; ct = params{21}; It = params{22};
ma = params{23}; ca = params{24}; Ia = params{25};
m0 = mb + mt + ma;

dt_sim = 0.002;
N_sim = 20000;
do_manip = true;
%% inertia
shape = [0 0 0; % <- x | y
         0 2 0; %      v
         0 0 0];
shape_mass = [mq, mq, mq;
              mq, m0, mq;
              mq, mq, mq];
shape_idx = shape;
[core_row, core_col] = find(shape == 2);
[AMs_rows, AMs_cols] = find(shape ~= 0);

num_AMs = length(AMs_cols);
AM_mass = 0; % mass of shape
Ir = zeros(3, 3); % inertia w.r.t. its com
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
    AM_com = AM_com + r_0j{j} * mass_ams(j)/AM_mass;
    shape_idx(AMs_rows(j), AMs_cols(j)) = j;
end
core_idx = shape_idx(core_row, core_col);

% compute Ir : tilting part
for j = 1:length(AMs_cols)
    r_cj{j} = -[(core_col - AMs_cols(j)) *-d ; (core_row - AMs_rows(j)) *d ;0] - AM_com;% p_com to j
    if j == core_idx
        Ij = It;
    else
        Ij = Iq;
    end
    I_cj{j} = Ij + mass_ams(j) * (r_cj{j}' * r_cj{j} * eye(3,3) - r_cj{j} * r_cj{j}'); % TODO
    Ir = Ir + I_cj{j};
end

if do_manip && norm(AM_com) < 1e-6
    fprintf("\nSymmetric shape : \n")
    disp(shape)
    fprintf("Do manipultion\n")
else
    fprintf("\nNot Symmetric shape : \n")
    disp(shape)
    fprintf("Do not manipultion\n")
end
%%
wn = 0.5; damp = 1.2; % 0.5, 1.2
kp_M = wn^2; 
kv_M = 2 * damp *sqrt(kp_M);

wn = 1.0; damp = 1.2; % 1, 1.2
kp_z = wn^2; 
kv_z = 2 * damp *sqrt(kp_z);

kp_M = diag([kp_M, kp_M, kp_z]);
kv_M = diag([kv_M, kv_M, kv_z]);

kw_I = 20; % 10 or 20
k_psi = 2; % 2 

%regulizer
k_lambda = 1.0;
k_delta = 1.0;

% servo moter
kp_servo = 0.01; % 0.01
kd_servo = 0.03; % 0.03
damp_servo = 0.1; % 0.1

% damped pseudoinverse
damped = 0.3; % 0.3
k_R = 1.0; % 1.0

epsilon = kv_M(1, 1) * 0.3; % kv_M(1, 1) * 0.3
alpha = 20; % 15 or 20
gamma = alpha * 3.0; % alpha * 3.0 
beta = diag([1 1 1]) * mq * norm(gravity) * 0.5; % 1 1 1 * mq * norm(gravity) * 0.5
N_sim_tmp = 10000;

mass_uncertainty = 1.10; 
inertia_uncertainty = 0.90;
thrust_limit = thrust_limit * 3;

%disturbance
sigma = 0.0; mean = 3.0; max_val = 30.0;
disturb_sim = zeros(N_sim, 6);
disturb = mean * [0; 0; 0; 0.5; -0.7; -1.0];
% X, w_estimation error
X_error = zeros(3, num_AMs);
w_error = zeros(3, num_AMs);
sigma_X = 3 / 100; max_X = 0.1;
sigma_w = 3 / 100; max_w = 0.1; 
delay = 0.01 / dt_sim; % position control delay dt_sim * delay
wz_delay = 0.1 / dt_sim;
delay_regul = 0.1 / dt_sim; % lambda/delta regularizer delay dt_sim * delay

rng('shuffle')z

if do_manip && norm(AM_com) < 1e-6
    X_hover = [1; 2; 3] * 1e-1;
    rpy_hover = [10, -20, 30] / 180 * pi;
    [X_des, Xd_des, Xdd_des, R_e_des, w_e_des] = get_traj_hover_manip(X_hover, rpy_hover, N_sim, dt_sim);

    % helix
    radius = 0.3;  v_z = 0.05;
    omega = 2 * pi * 0.1; 
    rpyd  = [0.00; -0.1; 0.1] * 2 * pi; 
    X_hover = [0; 0; 0.5]; rpy_hover = [10, -20, 30] / 180 * pi; 
    %[X_des, Xd_des, Xdd_des, R_e_des, w_e_des] = get_traj_helix_manip(radius, omega, v_z, rpyd, X_hover, rpy_hover, N_sim, dt_sim);
else
    X_hover = [1; 2; 3] * 1e-1; yaw_hover = 0 / 180 *pi; 
    [X_des, Xd_des, Xdd_des, yaw_des, yawd_des, yawdd_des] = get_traj_hover(X_hover, yaw_hover, N_sim, dt_sim);
    radius = 0.3;  v_z = 0.05;
    omega     = 2 * pi * 0.1; 
    omega_yaw = 2 * pi * 0.1; 
    X_hover = [0; 0; 0.5]; yaw_hover = 30 / 180 *pi; 
    %[X_des, Xd_des, Xdd_des, yaw_des, yawd_des, yawdd_des] = get_traj_helix(radius, omega, omega_yaw, v_z, X_hover, yaw_hover, N_sim, dt_sim);
end
%%
B = [r r -r -r;  -r r r -r; mu -mu mu -mu; 1 1 1 1];

e_1 = [1; 0; 0];
e_2 = [0; 1; 0];
e_3 = [0; 0; 1];
%%
X = [0; 0; 0]; Xd = [0; 0; 0]; Xdd = [0; 0; 0];
w = [0; 0; 0]; wd = [0; 0; 0];
theta = [0; 0]; thetad = [0; 0]; thetadd = [0; 0];
thetad_ref = [0; 0];
theta_ref  = [0; 0];
R = eye(3, 3);
Rd = R * S(w);

force_prev = mass_ams * norm(gravity) * 1.0;
delta_hat = zeros(3, num_AMs);
w_des_prev = zeros(3, num_AMs);

X_hist = []; Xd_hist = []; w_hist = []; wd_hist = []; R_hist = cell(N_sim, 1); 
w_des_hist = []; wd_des_hist = []; thrusts_hist = []; F_hat_hist = [];
e_p_hist = []; e_d_hist = []; e_w_hist = []; e_psi_hist = [];
nu_e_hist = []; delta_hat_hist = []; delta_tilde_hist = [];
force_per_M_hist = []; delta_hat_x_per_M_hist = [];
e_theta_hist = []; e_thetad_hist = [];
e_R_hist = [];
times = [];

for i = 1:N_sim_tmp
    % Semi-decentralized control
    tau_tot = zeros(3, 1);
    force_tot = 0;
    thrusts = [];
    force_per_M = [];
    delta_hat_x_per_M = [];
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
        
        if mod(i, delay) == 1
        e_p  = (1 - X_error(:, j)).*X - X_des(:, i); % could be computed using X_com = X_j - R' * r_cj
        e_pd = (1 - X_error(:, j)).*Xd - Xd_des(:, i);
        end
        e_pdd_hat = gravity + force_prev(j) / mj * R * e_3 - Xdd_des(:, i);
        %e_pdd_hat = Xdd - Xdd_des(:, i);
        
        if i < N_sim
            Xddd_des = (Xdd_des(:, i+1) - Xdd_des(:, i)) / dt_sim;
        end
        
        % regularizer
        if mod(i, delay_regul) == 1 && i > 1            
            delta_hat(:, j) = delta_hat(:, j) - k_delta * (delta_hat(:, j) - mj / AM_mass * delta_hat_sum);
            force_prev(j) = force_prev(j) - k_lambda * (force_prev(j) - mj / AM_mass * force_sum);
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
        if do_manip && norm(AM_com) < 1e-6 % TODO: delay or error
            Ij = I_cj{j};
            kw_j = kw_I * norm(Ij);
    
            J = [0 0 cos(theta(1)); 0 1 0; 1 0 -sin(theta(1))];
            J_dagger = J' * inv(J * J'+ damped^2 * eye(3,3));
            
            R_e = R * Ry(theta(1)) * Rx(theta(2));
            e_R = 0.5 * vee(R_e_des{i}' * R_e - R_e' * R_e_des{i});
            w_e = - k_R * e_R + R_e' * R_e_des{i} * w_e_des(:, i);
            sol = J_dagger * (Rx(theta(1)) * Rx(theta(2)) * w_e - [w_xj_des; w_yj_des; 0]);
            w_zj_des = sol(1);
            
            if mod(i, wz_delay) == 1
                w_zj_des_delayed = w_zj_des;
            end

            if j == core_idx
                w_des = [w_xj_des; w_yj_des; w_zj_des]; 
                if i > 1
                    wd_des = w_des - w_des_prev(:, j);
                else
                    wd_des = [0; 0; 0];
                end
                w_des_prev(:, j) = w_des;
        
                wj = (1 - w_error(:, j)).*w;
                e_w = wj - w_des;

                thetad_ref = sol(2:3);

                [G, C] = get_GC(Ir, mb, cb, Ib, ma, ca, Ia, wj, theta, thetad);
                tau_j = G(1:3, :) * [wd_des; thetadd] + C(1:3, :) * [wj; thetad] - kw_j * e_w;
 
                theta_ref = theta_ref + thetad_ref * dt_sim;
                e_theta = theta - theta_ref;
                e_thetad = thetad - thetad_ref;
                tau_theta = - kp_servo * e_theta - kd_servo * e_thetad;
            else
                w_des = [w_xj_des; w_yj_des; w_zj_des_delayed]; 
                if i > 1
                    wd_des = w_des - w_des_prev(:, j);
                else
                    wd_des = [0; 0; 0];
                end
                w_des_prev(:, j) = w_des;
        
                wj = (1 - w_error(:, j)).*w;
                e_w = wj - w_des;

                tau_j = S(wj) * Ij * wj + Ij * wd_des - kw_j * e_w;
            end
        else % do not manipulation
            psi = atan2(R(2,1), R(1,1));
            e_psi = psi - yaw_des(i);
            e_psi = wrapToPi(e_psi);
            psid_des = - k_psi * e_psi + yawd_des(i);
            w_zj_des = (psid_des - R(3,1) * w_xj_des - R(3,2) * w_yj_des ) / R(3,3);
            w_des = [w_xj_des; w_yj_des; w_zj_des]; 
            if i > 1
                wd_des = w_des - w_des_prev(:, j);
            else
                wd_des = [0; 0; 0];
            end
            w_des_prev(:, j) = w_des;
    
            wj = (1 - w_error(:, j)).*w;
            e_w = wj - w_des;
            Ij = I_cj{j};
            kw_j = kw_I * norm(Ij);
    
            if j == core_idx
                [G, C] = get_GC(Ir, mb, cb, Ib, ma, ca, Ia, wj, theta, thetad);
                tau_j = G(1:3, :) * [wd_des; thetadd] + C(1:3, :) * [wj; thetad] - kw_j * e_w;
                
                thetad_ref = [0; 0];
                theta_ref = theta_ref + thetad_ref * dt_sim;
                e_theta = theta - theta_ref;
                e_thetad = thetad - thetad_ref;
                tau_theta = - kp_servo * e_theta - kd_servo * e_thetad;
            else
                tau_j = S(wj) * Ij * wj + Ij * wd_des - kw_j * e_w;
            end
        end

        lambda_j = inv(B) * [tau_j; force_j];
        lambda_j = min(thrust_limit, max(-thrust_limit, lambda_j));
        
        tau_j = B(1:3,:) * lambda_j;
        force_j = B(4, :) * lambda_j;
        
        tau_tot = tau_tot + tau_j + S(r_cj{j}) * [0; 0; force_j];
        force_tot = force_tot + force_j;
        
        thrusts = [thrusts; lambda_j];
        force_per_M = [force_per_M; force_j / mj];
        delta_hat_x_per_M = [delta_hat_x_per_M; delta_hat(3, j) / mj];
    end
    delta_hat_sum = sum(delta_hat, 2);
    force_sum = force_tot;

    % generate disturbance
    disturb_dot = randn(6, 1) * sigma;
    disturb = disturb + disturb_dot * dt_sim;
    disturb = min(max(disturb, - max_val), max_val);
    disturb_sim(i,:) = disturb;

    % Dynamics
    tau_theta = tau_theta - damp_servo * thetad;
    [Xdd, wd, thetadd, Rd] = FD(AM_mass, Ir, mb, cb, Ib, ma, ca, Ia, R, w, theta, thetad,...
        force_tot, tau_tot, tau_theta, disturb, gravity, mass_uncertainty, inertia_uncertainty);
    
    Xd = Xd + Xdd * dt_sim;
    X = X + Xd * dt_sim + 0.5 * Xdd * dt_sim^2 ;
    w = w + wd * dt_sim;
    thetad = thetad + thetadd * dt_sim;
    theta = theta + thetad * dt_sim + 0.5 * thetadd * dt_sim^2 ;
    R = R + Rd * dt_sim;
    [U_, ~, V_] = svd(R);
    R  = U_ * V_';
    
    %thetad = thetad_ref;
    %theta = theta + thetad * dt_sim;

    % Log
    X_hist  = [X_hist, X];
    Xd_hist = [Xd_hist, Xd];
    w_hist  = [w_hist, w];
    w_des_hist = [w_des_hist, w_des];
    wd_hist = [wd_hist, wd];
    wd_des_hist = [wd_des_hist, wd_des]; 
    R_hist{i} = R;
    thrusts_hist = [thrusts_hist, thrusts];
    force_per_M_hist = [force_per_M_hist, force_per_M];
    delta_hat_x_per_M_hist = [delta_hat_x_per_M_hist, delta_hat_x_per_M];

    e_p = X - X_des(:, i); 
    e_pd = Xd - Xd_des(:, i);
    e_w = w - w_des;
   
    e_p_hist = [e_p_hist, e_p]; 
    e_d_hist = [e_d_hist, e_pd]; 
    e_w_hist = [e_w_hist, e_w];
    nu_e_hist = [nu_e_hist, nu_ej];
    delta_hat_hist = [delta_hat_hist, delta_hat(:, j)];
    delta_tilde_hist = [delta_tilde_hist, mj / AM_mass * disturb(4:6) - delta_hat(:, j)];

    e_theta_hist = [e_theta_hist, e_theta];
    e_thetad_hist = [e_thetad_hist, e_thetad]; 
    
    if do_manip && norm(AM_com) < 1e-6
        e_R_hist = [e_R_hist, e_R];
    else 
        e_psi_hist = [e_psi_hist, e_psi];
    end

    times = [times; i * dt_sim];
end
%% State plot
figure('Position',[50 350 700 600]);
colors = lines(6);

subplot(2,2,1)
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
subplot(2,2,2)
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
subplot(2,2,3)
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
if ~(do_manip && norm(AM_com) < 1e-6)
subplot(2,2,4)
hold on
plot(times, e_psi_hist, 'LineWidth', 1.0);
legend({'$e_\psi$'},'Interpreter', 'latex','FontSize', 12);
title('$e_\psi$', 'Interpreter', 'latex','FontSize', 14)
grid on
end
%% Error 
figure('Position',[750 350 600 600]);
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
%% Assumption: force/disturbance balance
figure('Position',[1350 550 500 400]);

% 7. thrusts
subplot(3,1,1)
hold on
plot(times, thrusts_hist, 'LineWidth', 1.0);
%ylim([ -thrust_limit, thrust_limit] )
title('Thrusts', 'Interpreter', 'latex','FontSize', 14)
grid on


legend_entries = arrayfun(@(x) sprintf('%d', x), 1:num_AMs, 'UniformOutput', false);
subplot(3,1,2)
plot(times, force_per_M_hist)
title('$\frac{\lambda_{p,i}}{M_i}$', 'Interpreter', 'latex','FontSize', 14)
legend(legend_entries)
grid on

subplot(3,1,3)
plot(times, delta_hat_x_per_M_hist)
title('$\frac{\Delta_{p,i}}{M_i}$', 'Interpreter', 'latex','FontSize', 14)
legend(legend_entries)
grid on
%% End effector Error 
figure('Position',[1350 50 500 450]);
colors = lines(6);

% e_R
if do_manip && norm(AM_com) < 1e-6
subplot(3,1,1)
hold on
for j = 1:3
    plot(times, e_R_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_R$', 'Interpreter', 'latex','FontSize', 14)
grid on
end

subplot(3,1,2)
hold on
for j = 1:2
    plot(times, e_theta_hist(j, :), 'Color', colors(3-j,:), 'LineWidth', 1.0);
end
legend({'$roll$','$pitch$'},'Interpreter','latex','FontSize', 12);
title('$e_\theta$', 'Interpreter', 'latex','FontSize', 14)
grid on

subplot(3,1,3)
hold on
for j = 1:2
    plot(times, e_thetad_hist(j, :), 'Color', colors(3-j,:), 'LineWidth', 1.0);
end
legend({'$roll$','$pitch$'},'Interpreter','latex','FontSize', 12);
title('$e_{\dot\theta}$', 'Interpreter', 'latex','FontSize', 14)
grid on


%%
function [Xdd, wd, thetadd, Rd] = FD(mass, Ir, mb, cb, Ib, ma, ca, Ia, R, w, theta, thetad, ...
            force, tau_w, tau_theta, disturb, gravity, mass_uncertainty, inertia_uncertainty)
    e_1 = [1; 0; 0];
    e_2 = [0; 1; 0];
    e_3 = [0; 0; 1];
    % model uncertainty
    mass = mass * mass_uncertainty;
    mb = mb * mass_uncertainty;
    ma = ma * mass_uncertainty;

    Ir = Ir * inertia_uncertainty;
    Ib = Ib * inertia_uncertainty;
    Ia = Ia * inertia_uncertainty;
    
    % position
    Xdd = (force * R * e_3 - mass * norm(gravity) * e_3 + disturb(4:6)) / mass;
    
    % rotaion
    J_b_1  = Ry(-theta(1));
    J_b_2  = [e_2, zeros(3,1)];
    r_b = - thetad(1) * e_2;
    J_b_1d = J_b_1 * S(r_b);

    J_a_1  = Rx(-theta(2)) * Ry(-theta(1));
    J_a_2  = [Rx(-theta(2)) * e_2, e_1];
    r_a = - thetad(2) * Ry(theta(1)) * e_1 - thetad(1) * e_2;
    J_a_1d = J_a_1 * S(r_a);
    J_a_2d = [Rx(-theta(2)) * S(-thetad(2) * e_1) * e_2, zeros(3,1)];

    I_b_bar = J_b_1' * Ib * J_b_1;
    I_b_tilde = J_b_1' * Ib * J_b_2;
    I_b_hat = J_b_2' * Ib * J_b_2;

    I_a_bar = J_a_1' * Ia * J_a_1;
    I_a_tilde = J_a_1' * Ia * J_a_2;
    I_a_hat = J_a_2' * Ia * J_a_2;
    
    I_c_bar = J_b_1' * (- ma * S2(ca) - mb * S2(cb)) * J_b_1;
    I_c_bard = S(-r_b) * I_c_bar + I_c_bar * S(r_b);
    I_c_hat = J_b_2' * (- ma * S2(ca) - mb * S2(cb)) * J_b_2;
    
    M_w = Ir + I_b_bar + I_a_bar + I_c_bar;
    M_w_theta = I_b_tilde + I_a_tilde;
    M_theta = I_b_hat + I_a_hat + I_c_hat;

    C_w = S(w) * Ir + S(w - r_b) * I_b_bar + S(w - r_b) * I_a_bar + ...
          I_b_bar * S(r_b) + I_a_bar * S(r_b) + I_c_bard + S(w) * I_c_bar;
    C_w_theta = S(w - r_b) * I_b_tilde + S(w - r_b) * I_a_tilde + J_b_1'* Ia * J_a_2d;
    C_theta_w = J_a_2d'* Ia * J_a_1 + I_a_tilde' * S(r_a) + I_b_tilde' * S(r_b);
    C_theta = J_a_2d'* Ia * J_a_2 + J_a_2'* Ia * J_a_2d;

    G = [M_w, M_w_theta; M_w_theta' M_theta];
    C = [C_w, C_w_theta; C_theta_w, C_theta];
    V = inv(G) * ([tau_w; tau_theta] - C * [w; thetad] + [disturb(1:3); 0; 0]);
    wd = V(1:3);
    thetadd = V(4:5);
    Rd = R * S(w);
end
function [G, C] = get_GC(Ir, mb, cb, Ib, ma, ca, Ia, w, theta, thetad)
    e_1 = [1; 0; 0];
    e_2 = [0; 1; 0];
    e_3 = [0; 0; 1];
 
    % rotaion
    J_b_1  = Ry(-theta(1));
    J_b_2  = [e_2, zeros(3,1)];
    r_b = - thetad(1) * e_2;
    J_b_1d = J_b_1 * S(r_b);

    J_a_1  = Rx(-theta(2)) * Ry(-theta(1));
    J_a_2  = [Rx(-theta(2)) * e_2, e_1];
    r_a = - thetad(2) * Ry(theta(1)) * e_1 - thetad(1) * e_2;
    J_a_1d = J_a_1 * S(r_a);
    J_a_2d = [Rx(-theta(2)) * S(-thetad(2) * e_1) * e_2, zeros(3,1)];

    I_b_bar = J_b_1' * Ib * J_b_1;
    I_b_tilde = J_b_1' * Ib * J_b_2;
    I_b_hat = J_b_2' * Ib * J_b_2;

    I_a_bar = J_a_1' * Ia * J_a_1;
    I_a_tilde = J_a_1' * Ia * J_a_2;
    I_a_hat = J_a_2' * Ia * J_a_2;
    
    I_c_bar = J_b_1' * (- ma * S2(ca) - mb * S2(cb)) * J_b_1;
    I_c_bard = S(-r_b) * I_c_bar + I_c_bar * S(r_b);
    I_c_hat = J_b_2' * (- ma * S2(ca) - mb * S2(cb)) * J_b_2;
    
    M_w = Ir + I_b_bar + I_a_bar + I_c_bar;
    M_w_theta = I_b_tilde + I_a_tilde;
    M_theta = I_b_hat + I_a_hat + I_c_hat;

    C_w = S(w) * Ir + S(w - r_b) * I_b_bar + S(w - r_b) * I_a_bar + ...
          I_b_bar * S(r_b) + I_a_bar * S(r_b) + I_c_bard + S(w) * I_c_bar;
    C_w_theta = S(w - r_b) * I_b_tilde + S(w - r_b) * I_a_tilde + J_b_1'* Ia * J_a_2d;
    C_theta_w = J_a_2d'* Ia * J_a_1 + I_a_tilde' * S(r_a) + I_b_tilde' * S(r_b);
    C_theta = J_a_2d'* Ia * J_a_2 + J_a_2'* Ia * J_a_2d;

    G = [M_w, M_w_theta; M_w_theta' M_theta];
    C = [C_w, C_w_theta; C_theta_w, C_theta];
end

function out = Rx(q)
    out = [1, 0, 0;
          0, cos(q), -sin(q);
          0, sin(q), cos(q)];
end
function out = Ry(q)
    out = [cos(q), 0, sin(q);
           0, 1, 0;
           -sin(q), 0, cos(q)];
end