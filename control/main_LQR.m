addpath("../../casadi-3.6.7-windows64-matlab2018b" , "dynamics", "casadi_functions", "functions", "../params" )
clear; load("../3D_ver2/data/test.mat")
close all;
%%
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

N_horizon = 10; do_print = 0; do_plot = 1;
[min_eigval_arr, max_eigval_arr, xT_W_r_x_arr] = check_rechability_gramian(A_nom_cell, B_nom_cell, N_horizon, dt_sim, N_sim, n, do_print, do_plot);

Q  = diag([1,1,1,1,1,1,1,1]);
Qf = Q;
R  = eye(nu) * 1e-1;
K_arr = compute_LQR_gains(A_nom_cell, B_nom_cell, Q, Qf, R, dt_sim);
%% Simulation
delta_inertia_nom = 1.0; delta_k_nom = 1.0;
delta_inertia = 1.0; delta_k = 1.0;
sigma = 0.0; mean = 0.00; max_val = 0.5;

x_sim = zeros(N_sim + 1, nx); x_sim(1,:) = x_d_interp(1,:);
u_sim = zeros(N_sim, nu);

disturb_sim = zeros(N_sim, n);
disturb = mean * ones(n,1);

rng('shuffle')
for i = 1:N_sim
    disturb_dot = randn(n, 1) * sigma;
    disturb = disturb + disturb_dot * dt_sim;
    disturb = min(max(disturb, - max_val), max_val);
    disturb_sim(i,:) = disturb;

    u_fb = u_d_interp(i,:)' - K_arr{i}*(x_sim(i,:)' - x_d_interp(i,:)');

    x_dot_act = full( x_dot_func(x_sim(i,:), u_fb, delta_inertia, delta_k, disturb) );

    u_sim(i, :) = u_fb';

    x_sim(i+1,:) = x_sim(i,:) + x_dot_act' * dt_sim;
end


%% plot
plot_simulation_results(t_sim, x_sim, x_d_interp, u_sim, u_d_interp, disturb_sim, shape, dh, gravity, ...
    n, nx, nu, N_sim, dt_sim, delta_inertia, delta_k, sigma, mean, max_val)

% plot video
[AM_com, AM_mass, AM_inertia] = get_inertia(shape, m0, I0, d);
mass = {mass_door(1), mass_door(2), mass_door(3), AM_mass};
inertia{4} = AM_inertia;
r_i_ci{4} = [AM_com(1); 0; AM_com(2)];

slow_factor = 1; force_scale = 0.2; do_view = 0; q = [0; 0; 0; 0];
robot = generate_door_ver2(n, dh, r_i_ci, d, gravity, shape, mass, inertia, do_view, q);

%save_plot_tree(robot, dh, params, x_sim, u_sim, dt_sim, N_sim, slow_factor, force_scale, shape)
plot_tree(robot, dh, params, x_sim, u_sim, dt_sim, N_sim, slow_factor, force_scale, shape)