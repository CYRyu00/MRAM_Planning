addpath("../../casadi-3.6.7-windows64-matlab2018b" , "dynamics", "casadi_functions", "functions", "../params" )
clear; 
load("../planning_continous/data/Q2_1e0_1110/10sec/6_1_0.mat")
close all;
%% CONTROL GAIN
alpha = 1e-2; beta = 1e-4; % 1e-2 / 1e-4
b = 2e-1; k = 3e-1; % 2e-1/ 3e -1
do_video = true;
%% Parsing and Interpolation 
shape = zeros([K,L]);
x_d = x_opt; u_d = zeros([N, 4 * num_AMs]);
nx = 2 * n; nu = 4 * num_AMs;
[AMs_rows, AMs_cols] = find(rho_opt >= 0.9);
for i = 1:length(AMs_rows)
    shape(AMs_rows(i), AMs_cols(i)) = 1;
    u_d(:, 4*i-3:4*i) = u_opt(:, 4*((AMs_rows(i)-1)*L + AMs_cols(i)) - 1 : 4 * ((AMs_rows(i)-1)*L + AMs_cols(i)) + 2);
end
shape(core(1), core(2)) = 2;
fprintf('number of AMs : %d, Shape : \n', num_AMs)
disp(shape)

params = define_params();
[x_dot_func, A_func, B_func] = define_dynamics(shape, num_AMs, params);

dt_sim = 0.01;
N_sim = N * dt / dt_sim;
t_plan = linspace(0, dt * N, N + 1);
t_sim = linspace(0, dt_sim * N_sim, N_sim + 1);

do_plot = 0;
[x_d_interp, u_d_interp] = interpolate_traj(x_d, u_d, t_plan, t_sim, do_plot);

delta_inertia = 1.0; delta_k = 1.0; disturb = [0; 0; 0; 0];
[A_nom_cell, B_nom_cell] = check_ctrb(x_d_interp, u_d_interp, N_sim, A_func, B_func, delta_inertia, delta_k, disturb);

N_horizon = 10; do_print = 0; do_plot = 0;
[min_eigval_arr, max_eigval_arr, xT_W_r_x_arr] = check_rechability_gramian(A_nom_cell, B_nom_cell, N_horizon, dt_sim, N_sim, n, do_print, do_plot);
%% Simulation
delta_inertia = 1.0; delta_k = 1.0;
sigma = 0.0; mean = 0.10; max_val = 0.0;

x_sim = zeros(N_sim + 1, nx); x_sim(1,:) = x_d_interp(1,:);
u_sim = zeros(N_sim, nu);

Y_err = zeros(N_sim, 3);
Yd_err = zeros(N_sim, 3);
w_err = zeros(N_sim, 3);
wd_err = zeros(N_sim, 3);
wd_history = zeros(N_sim, 3);
wd_des_history = zeros(N_sim, 3);

disturb_sim = zeros(N_sim, n);
disturb = mean * [1; 1; 1; 1];
rng('shuffle')

[core_row, core_col] = find(shape == 2);
[AMs_rows, AMs_cols] = find(shape ~= 0);
    
A = [r r -r -r;  -r r r -r;mu -mu mu -mu; 1 1 1 1];

int_wd_des = [0; 0; 0];
u_fb = [];

Y = cell(num_AMs, 1);
Yd = cell(num_AMs, 1);
Y_des = cell(num_AMs, 1);
Yd_des = cell(num_AMs, 1); 
Yd_des_prev = cell(num_AMs, 1);
for j = 1:length(AMs_rows)
    Yd_des_prev{j} = zeros(3, 1);
end

for i = 1:N_sim
    % generate disturbance
    disturb_dot = randn(n, 1) * sigma;
    disturb = disturb + disturb_dot * dt_sim;
    disturb = min(max(disturb, - max_val), max_val);
    disturb_sim(i,:) = disturb;
    
    % get simulated states
    [T_0i, w_cell, v_cell] = get_T_w_v(x_sim(i, 1:n)', x_sim(i, n+1:end)', dh, n);
    R_45 = [1 0 0;0 0 -1;0 1 0];
    R = T_0i{4}(1:3, 1:3) * R_45;
    X = T_0i{4}(1:3, 4);
    Xd = v_cell{4};
    w = R_45 * w_cell{4};
    
    % get desired states
    [T_0i_des, w_cell_des, v_cell_des] = get_T_w_v(x_d_interp(i, 1:n)', x_d_interp(i, n+1:end)', dh, n);
    R_des = T_0i_des{4}(1:3, 1:3) * R_45;
    X_des = T_0i_des{4}(1:3, 4);
    Xd_des = v_cell_des{4};
    w_des = R_45 * w_cell_des{4};
    
    % calculate decentralized control thrusts
    for j = 1:length(AMs_rows)
        rj = [(core_col - AMs_cols(j)) *-d ; (core_row - AMs_rows(j)) *d ; -r_i_ci{4}(2)]; % p_j,core
        dj = [0.00; 0; 0.03]; % [0.00; 0; 0.03] 
        Y{j} = X + R * (dj + rj);
        Yd{j} = Xd + R * S(w) * (dj + rj);

        Y_des{j} = X_des + R_des * (dj + rj);
        Yd_des{j} = Xd_des + R_des * S(w_des) * (dj + rj);
        Ydd_des = (Yd_des{j} - Yd_des_prev{j}) / dt_sim;
        Yd_des_prev{j} = Yd_des{j}; 

        uj = m0 * Ydd_des - b * (Yd{j} - Yd_des{j}) - k * (Y{j} - Y_des{j});

        nu_ = R' * (uj + m0 * norm(gravity) * [0; 0; 1]) - m0 * S2(w) *dj;
        wd_des = [-nu_(2)/m0/dj(3); nu_(1)/m0/dj(3); 0];
        lambda = nu_(3) + dj(1) / dj(3) * nu_(1);

        int_wd_des = int_wd_des + wd_des * dt_sim;
        tau_ = S(w) * (I0 * w) + I0*(wd_des - alpha *(w - int_wd_des)) - beta * w;

        u_fb(4 * j - 3 : 4 * j) = A\[tau_; lambda];
    end

    x_dot = full( x_dot_func(x_sim(i,:), u_fb, delta_inertia, delta_k, disturb) );
    
    u_sim(i, :) = u_fb';
    x_sim(i+1,:) = x_sim(i,:) + x_dot' * dt_sim;

    [wd_cell] = get_wd(x_sim(i, 1:n)', x_sim(i, n+1:end)', x_dot(n+1:end)', dh, n);
    wd = wd_cell{4};
    j = 1; %you could chose 1 ~ num_AMs
    Yd_err(i, :) = Yd{j} - Yd_des{j};
    Y_err(i, :) = Y{j} - Y_des{j};
    w_err(i, :) = w - int_wd_des;
    wd_err(i, :) = wd - wd_des;
    wd_history(i, :) = wd';
    wd_des_history(i, :) = wd_des';
end


%% plot
close all
plot_simulation_results(t_sim, x_sim, x_d_interp, u_sim, u_d_interp, disturb_sim, shape, dh, gravity, ...
    n, nx, nu, N_sim, dt_sim, delta_inertia, delta_k, sigma, mean, max_val)

colors = lines(3);
figure('Position',[600, 50, 800, 700])
times = (1:N_sim) * dt_sim;

% LaTeX-compatible subplot titles and legends
subplot(3,2,1)
plot(times, Y_err)
legend({'$Y_1$', '$Y_2$', '$Y_3$'}, 'Interpreter','latex', 'FontSize',10)
title('${Y}$ error', 'Interpreter','latex','FontSize', 14)

subplot(3,2,2)
plot(times, Yd_err)
legend({'$\dot{Y}_1$', '$\dot{Y}_2$', '$\dot{Y}_3$'}, 'Interpreter','latex', 'FontSize',10)
title('$\dot{{Y}}$ error', 'Interpreter','latex','FontSize', 14)

subplot(3,2,3)
plot(times, w_err)
legend({'$\omega_1$', '$\omega_2$', '$\omega_3$'}, 'Interpreter','latex', 'FontSize',10)
title('$\mathbf{\omega}$ error', 'Interpreter','latex','FontSize', 14)

subplot(3,2,4)
plot(times, wd_err)
legend({'$\dot{\omega}_1$', '$\dot{\omega}_2$', '$\dot{\omega}_3$'}, 'Interpreter','latex', 'FontSize',10)
title('$\dot{\mathbf{\omega}}$ error', 'Interpreter','latex','FontSize', 14)

subplot(3,2,5)
hold on
colors = lines(3); % Ensure colors is defined
for i = 1:3
    plot(times, wd_history(:,i), 'Color', colors(i,:), 'LineWidth', 1.5)     
end
legend({'$\dot{\omega}_1$', '$\dot{\omega}_2$', '$\dot{\omega}_3$'}, ...
       'Interpreter','latex', 'FontSize', 10);
title('$\dot{\mathbf{\omega}}$', 'Interpreter', 'latex', 'FontSize', 14)

subplot(3,2,6)
hold on
colors = lines(3); % Ensure colors is defined
for i = 1:3
    plot(times, wd_history(:,i), 'Color', colors(i,:), 'LineWidth', 1.5)           
    plot(times, wd_des_history(:,i), '--', 'Color', colors(i,:), 'LineWidth', 1)
end
legend({'$\dot{\omega}_1$', '$\dot{\omega}_1^{\mathrm{des}}$', ...
        '$\dot{\omega}_2$', '$\dot{\omega}_2^{\mathrm{des}}$', ...
        '$\dot{\omega}_3$', '$\dot{\omega}_3^{\mathrm{des}}$'}, ...
       'Interpreter', 'latex', 'FontSize', 10);
title('$\dot{\mathbf{\omega}}$ vs $\dot{\mathbf{\omega}}_{\mathrm{des}}$', 'Interpreter','latex','FontSize', 14)


if do_video
    [AM_com, AM_mass, AM_inertia] = get_inertia(shape, m0, I0, d);
    mass = {mass_door(1), mass_door(2), mass_door(3), AM_mass};
    inertia{4} = AM_inertia;
    r_i_ci{4} = [AM_com(1); r_i_ci{4}(2); AM_com(2)];
    
    slow_factor = 1; force_scale = 0.2; do_view = 0; q = [0; 0; 0; 0];
    robot = generate_door_ver2(n, dh, r_i_ci, d, gravity, shape, mass, inertia, do_view, q);
    
    %save_plot_tree(robot, dh, params, x_sim, u_sim, dt_sim, N_sim, slow_factor, force_scale, shape)
    plot_tree(robot, dh, params, x_sim, u_sim, dt_sim, N_sim, slow_factor, force_scale, shape)
end