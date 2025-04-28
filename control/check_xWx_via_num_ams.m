addpath("../../casadi-3.6.7-windows64-matlab2018b" , "dynamics", "../params")
num_AMs_arr = 1:5;
min_eigval_arr = [];
xT_W_r_x_cell = cell(1, length(num_AMs_arr));

for num_AMs = num_AMs_arr
    filename = sprintf('../3D_ver2/data/result_9_5/ref_1/%d_1_0.mat', num_AMs);
    load(filename)
    shape = zeros([K, L]);
    shape_idx = zeros([K, L]);
    x_d = x_opt;
    u_d = zeros([N, 4 * num_AMs]);
    [AMs_rows, AMs_cols] = find(rho_opt >= 0.9);
    for i = 1:length(AMs_rows)
        shape(AMs_rows(i), AMs_cols(i)) = 1;
        shape_idx(AMs_rows(i), AMs_cols(i)) = i;
        u_d(:, 4 * i - 3 : 4 * i) = u_opt(:, 4 * ((AMs_rows(i)-1) * L + AMs_cols(i)) - 1 : 4 * ((AMs_rows(i) - 1) * L + AMs_cols(i)) + 2 );
    end
    shape(core(1), core(2)) = 2;
    fprintf('num AMs: %d, Shape : \n', num_AMs)
    disp(shape)

    params = define_params;
    [x_dot_func, A_func, B_func] = define_dynamics(shape, num_AMs, params);

    do_plot = 0;
    [x_d_interp, u_d_interp] = interpolate_traj(x_d, u_d, t_plan, t_sim, do_plot);

    delta_inertia = 1.0; delta_k = 1.0; disturb = [0; 0; 0; 0];
    [A_nom_cell, B_nom_cell] = check_ctrb(x_d_interp, u_d_interp, N_sim, A_func, B_func, delta_inertia, delta_k, disturb);

    N_horizon = 10; do_print = 0; do_plot = 0;
    [min_eigvalues, max_eigvalues, xT_W_r_x_arr ] = check_rechability_gramian(A_nom_cell, B_nom_cell, N_horizon, dt_sim, N_sim, n, do_print, do_plot);
    min_eigval_arr = [min_eigval_arr; min_eigvalues'];
    xT_W_r_x_cell{num_AMs} = xT_W_r_x_arr;
end
%% plot
close all
t_tmp = (1:N_sim - N_horizon + 1) * dt_sim;
legend_entries = arrayfun(@(k) ['$' num2str(k) '$'], num_AMs_arr, 'UniformOutput', false);

figure
plot(t_tmp, log(min_eigval_arr) / log(10))
legend(legend_entries, 'FontSize', 10, 'Interpreter', 'latex');
title('minimum eigenvalue for each number of AM' ...
    , 'FontSize', 14, 'Interpreter', 'latex');
xlabel('time[sec]', 'FontSize', 14, 'Interpreter', 'latex')
ylabel( '$\log_{10}{(\min{\lambda})}$','FontSize', 14,'Interpreter', 'latex')
grid on

figure
hold on
for num_AMs = num_AMs_arr
    plot(t_tmp, log(xT_W_r_x_cell{num_AMs}(:, 1)) / log(10))
end
legend(legend_entries, 'FontSize', 10, 'Interpreter', 'latex');
title({'$e_1^T W_r e_1$ for each number of AM'}, ...
    'Interpreter','latex','FontSize',14);
xlabel('time[sec]', 'FontSize', 14, 'Interpreter', 'latex')
ylabel( '$\log_{10}{(e_1^T W_r e_1)}$','FontSize', 14,'Interpreter', 'latex')
grid on