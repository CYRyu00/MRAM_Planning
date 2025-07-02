addpath("../../casadi-3.6.7-windows64-matlab2018b" , "dynamics", "casadi_functions", "functions", "../params" )
clear; close all;
load("../planning_continous_tilted/data/Q2_1e0_1110/2_1_0.mat")

m_duo = 2 * m0;
D = [0.5 * d; 0; 0];
I_duo = 2 * I0 + 2 * m0 * (D' * D * eye(3,3) - D * D');
%% CONTROL GAIN
wn = 2.0; damp = 1.2; % 3 / 1.2
k_p_x = m_duo * wn^2;
k_d_x = 2 * damp * sqrt(m_duo * k_p_x);  

wn = 2.0; damp = 1.2; % 3 / 1.2
k_p_z = m_duo * wn^2;
k_d_z = 2 * damp * sqrt(m_duo * k_p_z);  

k_p = diag([k_p_x, k_p_x, k_p_z]);
k_d = diag([k_d_x, k_d_x, k_d_z]);

wn = 1.0; damp = 1.2; % 1 / 1.2
k_R = m_duo * wn^2;
k_w = 2 * damp * sqrt(m_duo * k_R);  

do_video = true;
save_video = false;
N_sim_tmp = 1000;
dN = 10;
%% Parsing and Interpolation 
shape = zeros([K, L]);
x_d = x_opt; u_d = zeros([N, 8 * num_AMs]);
nx = 2 * n; nu = 8 * num_AMs;
[AMs_rows, AMs_cols] = find(rho_opt >= 0.9);
for i = 1:length(AMs_rows)
    shape(AMs_rows(i), AMs_cols(i)) = 1;
    u_d(:, 8*i-7:8*i) = u_opt(:, 8*((AMs_rows(i)-1)*L + AMs_cols(i)) - 7 : 8 * ((AMs_rows(i)-1)*L + AMs_cols(i)));
end
shape(core(1), core(2)) = 2;
fprintf('number of AMs : %d, Shape : \n', num_AMs)
disp(shape)

params = define_params();
x_dot_func = define_dynamics_duo(shape, num_AMs, params, theta);

dt_sim = 0.01;
N_sim = N * dt / dt_sim;
t_plan = linspace(0, dt * N, N + 1);
t_sim = linspace(0, dt_sim * N_sim, N_sim + 1);

do_plot = 0;
[x_d_interp, u_d_interp] = interpolate_traj(x_d, u_d, t_plan, t_sim, do_plot);
%% Simulation
delta_inertia = 1.0; delta_k = 1.0;
sigma = 0.0; mean = 0.00; max_val = 0.0;

x_sim = zeros(N_sim + 1, nx); x_sim(1,:) = x_d_interp(1,:);
u_sim = zeros(N_sim, nu);

disturb_sim = zeros(N_sim, n);
disturb = mean * [1; 1; 1; 1];
rng('shuffle')

[core_row, core_col] = find(shape == 2);
[AMs_rows, AMs_cols] = find(shape ~= 0);
    
s = sin(theta); c = cos(theta);
A_theta = [mu*s/sqrt(2) - r*c, mu*s/sqrt(2) - r*c, - mu*s/sqrt(2) + r*c, - mu*s/sqrt(2) + r*c;
     - mu*s/sqrt(2) + r*c, mu*s/sqrt(2) - r*c, mu*s/sqrt(2) - r*c, - mu*s/sqrt(2) + r*c;
     mu*c + sqrt(2)*r*s, - mu*c - sqrt(2)*r*s, mu*c + sqrt(2)*r*s, - mu*c - sqrt(2)*r*s;
     s/sqrt(2), - s/sqrt(2), - s/sqrt(2), s/sqrt(2);
     - s/sqrt(2), - s/sqrt(2), s/sqrt(2), s/sqrt(2);
     c, c, c, c];

p1 = [-0.5 * d; 0; 0]; p2 = [0.5 * d; 0; 0]; % p(:,3) = [D; 0; 0];
R1 = eye(3,3); R2 = [0 -1 0;1 0 0; 0 0 1]; % R{3} = eye(3,3);

A_duo = [Ad(R1, p1) * A_theta, Ad(R2, p2) * A_theta];
A_dagger = A_duo' * inv(A_duo * A_duo'); % least square


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

e_R = cell(num_AMs, 1);
e_w = cell(num_AMs, 1);
e_p = cell(num_AMs, 1);
e_d = cell(num_AMs, 1);

e_R_hist  = [];
e_w_hist  = [];
e_p_hist  = [];
e_d_hist = [];
times = [];

for i = 1:N_sim_tmp %N_sim
    % generate disturbance
    disturb_dot = randn(n, 1) * sigma;
    disturb = disturb + disturb_dot * dt_sim;
    disturb = min(max(disturb, - max_val), max_val);
    disturb_sim(i,:) = disturb;
    
    %TODO : check get T w v function 
    % get simulated states
    [T_0i, w_cell, v_cell] = get_T_w_v(x_sim(i, 1:n)', x_sim(i, n+1:end)', dh, n);
    R_45 = [1 0 0; 0 0 -1;0 1 0];
    R = T_0i{4}(1:3, 1:3) * R_45;
    X = T_0i{4}(1:3, 4);
    Xd = v_cell{4};
    w = R_45' * w_cell{4};
    
    % get desired states
    [T_0i_des, w_cell_des, v_cell_des] = get_T_w_v(x_d_interp(i, 1:n)', x_d_interp(i, n+1:end)', dh, n);
    R_des = T_0i_des{4}(1:3, 1:3) * R_45;
    X_des = T_0i_des{4}(1:3, 4);
    Xd_des = v_cell_des{4};
    w_des = R_45' * w_cell_des{4};
    if i == 1
        w_des_prev = w_des;
    end
    wd_des = (w_des - w_des_prev) / dt_sim;
    w_des_prev = w_des;
    
    % calculate decentralized control thrusts
    for j = 1:length(AMs_rows)
        rj = [(core_col - AMs_cols(j)) * -2 * d ; (core_row - AMs_rows(j)) *d ; -r_i_ci{4}(2)]; % p_j,core
        dj = [0.00; 0; 0.00]; % [0.00; 0; 0.03] 
        Y{j} = X + R * (dj + rj);
        Yd{j} = Xd + R * S(w) * (dj + rj);

        Y_des{j} = X_des + R_des * (dj + rj);
        Yd_des{j} = Xd_des + R_des * S(w_des) * (dj + rj);
        Ydd_des = (Yd_des{j} - Yd_des_prev{j}) / dt_sim;
        Yd_des_prev{j} = Yd_des{j}; 

        e_R{j} = 0.5 * vee(R_des' * R - R' * R_des);
        e_w{j} = w - R' * R_des * w_des;

        e_p{j} = Y{j} - Y_des{j};
        e_d{j} = Yd{j} - Yd_des{j};
        
        force_ = R' * ( m_duo * gravity + m_duo * Ydd_des - k_d * e_d{j} - k_p * e_p{j});
        tau_   = S(w) * I_duo * w - k_R * e_R{j} - k_w * e_w{j} ... % TODO e_I for integral
                 - I_duo * (S(w) * R_des * w_des - R_des * wd_des);
        u_fb(8 * j - 7 : 8 * j) = A_dagger * [tau_; force_];
    end
    %u_fb = u_d_interp(i, :);
    x_dot = full(x_dot_func(x_sim(i,:), u_fb, delta_inertia, delta_k, disturb));
    
    u_sim(i, :) = u_fb';
    x_sim(i+1,:) = x_sim(i,:) + x_dot' * dt_sim;

    [wd_cell] = get_wd(x_sim(i, 1:n)', x_sim(i, n+1:end)', x_dot(n+1:end)', dh, n);
    wd = wd_cell{4};

    j = 1; %you could chose 1 ~ num_AMs
    e_R_hist = [e_R_hist, e_R{j}];
    e_w_hist = [e_w_hist, e_w{j}];
    e_p_hist = [e_p_hist, e_p{j}];
    e_d_hist = [e_d_hist, e_d{j}];
    times = [times, i * dt_sim];
end


%% plot
close all
plot_simulation_results(t_sim, x_sim, x_d_interp, u_sim, u_d_interp, disturb_sim, shape, dh, gravity, ...
    n, nx, nu, N_sim, dt_sim, delta_inertia, delta_k, sigma, mean, max_val)

colors = lines(3);
figure('Position',[600, 50, 800, 700])

subplot(2,2,1)
plot(times, e_R_hist)
legend({'$x$', '$y$', '$z$'}, 'Interpreter','latex', 'FontSize',10)
title('$e_R$', 'Interpreter','latex','FontSize', 14)
grid on

subplot(2,2,2)
plot(times, e_w_hist)
legend({'$x$', '$y$', '$z$'}, 'Interpreter','latex', 'FontSize',10)
title('$e_w$', 'Interpreter','latex','FontSize', 14)
grid on

subplot(2,2,3)
plot(times, e_p_hist)
legend({'$x$', '$y$', '$z$'}, 'Interpreter','latex', 'FontSize',10)
title('$e_p$', 'Interpreter','latex','FontSize', 14)
grid on

subplot(2,2,4)
plot(times, e_d_hist)
legend({'$x$', '$y$', '$z$'}, 'Interpreter','latex', 'FontSize',10)
title('$e_{d}$', 'Interpreter','latex','FontSize', 14)
grid on

if do_video
    [AM_com, AM_mass, AM_inertia] = get_inertia_duo(shape, m0, I0, d);
    mass = {mass_door(1), mass_door(2), mass_door(3), AM_mass};
    inertia{4} = AM_inertia;
    r_i_ci{4} = [AM_com(1); r_i_ci{4}(2); AM_com(2)];
    
    slow_factor = 1.0; force_scale = 0.2; do_view = 0; q = [0; 0; 0; 0];
    robot = generate_door_duo(n, dh, r_i_ci, d, gravity, shape, mass, inertia, do_view, q);
    
    if save_video
        save_plot_tree(robot, dh, params, x_sim, u_sim, dt_sim, N_sim, slow_factor, force_scale, shape)
    else
        plot_tree_duo(robot, dh, params, theta, x_sim, u_sim, dt_sim, N_sim, slow_factor, force_scale, shape, dN)
    end
end