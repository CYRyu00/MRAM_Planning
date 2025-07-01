addpath("dynamics", "../params", "plot", "../../casadi-3.6.7-windows64-matlab2018b")
clear; close all;
%%
% Define Dynamic parameters
params = define_params();
m0 = params{1}; I0 = params{2}; mu = params{3}; r = params{4}; d = params{5};
thrust_limit= params{6}; kt = params{7}; c_1 = params{8}; c_2 = params{9}; mass_door = params{10};
handle_factor = params{11}; inertia = params{12}; r_i_ci = params{13}; n = params{14}; dh = params{15}; gravity = params{16};
theta = 15 / 180 * pi;

% NLP parameters
dt = 0.1; N = 100;

x_0 = [0; 0; 0; 0; 0; 0; 0; 0];
x_f = [pi/4; 0; 0; 0; 0; 0; 0; 0];

% Generate reference trajectory & initial guess of x
t0 = 1; % cosine smoothing 
t1 = 1; % hovering time
q_o_ref = generate_q_o_ref(x_0, x_f, N, dt, t0 ,t1);
x_interp = generate_x_interp(x_0, x_f, N, dt, t1);
%%
thrust_scale = 1;
tau_scale = 0;

u_max = thrust_limit * thrust_scale;
u_min = thrust_limit * (-thrust_scale);
tau_min = -0.2 * tau_scale ;
tau_max =  0.2 * tau_scale ;

max_iter = 2000;
eps = 0.25;
gamma = 0.3;

for num_AMs = 4:1:6
    for CASE = -1
        if CASE == 1 && num_AMs > 10
            continue
        end
        if CASE == -1
            K = 2 * num_AMs - 1; L = num_AMs; core = [num_AMs, 1];
            if num_AMs > 5
                K = 9; L = 5; core = [5, 1];
            end
            mkdir data/Q2_1e0_1110/10sec
            filename = sprintf('data/Q2_1e0_1110/10sec/%d_%d_%d.mat', num_AMs, thrust_scale, tau_scale);
        elseif CASE == 1
            K = 2 * num_AMs - 1; L = num_AMs; core = [num_AMs, 1];
            mkdir data/result_full/ref_2
            filename = sprintf('data/result_full/ref_2/%d_%d_%d.mat', num_AMs, thrust_scale, tau_scale);
        elseif CASE == 2
            K = 9; L = 5; core = [5, 1];
            mkdir data/result_9_5/tool
            filename = sprintf('data/result_9_5/tool/%d_%d_%d.mat', num_AMs, thrust_scale, tau_scale);
        elseif CASE == 3
            K = 11; L = 6; core = [6, 1];
            mkdir data/result_11_6/ref_2
            filename = sprintf('data/result_11_6/ref_2/%d_%d_%d.mat', num_AMs, thrust_scale, tau_scale);
        end
        fprintf("CASE = %d, [K, L] = [%d , %d]\n", CASE, K, L);
        fprintf("max_iter = %d, eps = %.2f, gamma = %.2f\n", max_iter, eps, gamma);
        fprintf("file name : "); disp(filename);

        nu = K * L * 8; zero_us = zeros(nu, 1);
        rho_init = ones(K, L) / K / L * (num_AMs - 1);
        rho_init(core(1), core(2)) = 1;
        X_init_guess = [reshape(rho_init, K * L, 1); reshape(x_interp', (N + 1) * 8, 1); repmat(zero_us, N, 1)];

        [Q1, Q2, Q3, R] = define_nlp_params_ver2(nu);

        %Solve NLP
        [rho_opt, x_opt, u_opt, optimal_value, exit_flag, processing_time, rho_history, exit_flag_history, time_history] ....
            =  solve_nlp(params, theta, num_AMs, K, L, core, q_o_ref, tau_min, tau_max, ...
                         u_min, u_max, x_0, x_f, X_init_guess, dt, N, t1, max_iter, eps, gamma, Q1, Q2, Q3, R);
        fprintf("num AMs: %d\n", num_AMs);
        fprintf("exit flag: %d \n", exit_flag);
        fprintf("optimal value: %f \n", optimal_value);
        fprintf("rho: \n"); disp(rho_opt)

        save(filename);
    end
end
%% plot
plot_graph_video