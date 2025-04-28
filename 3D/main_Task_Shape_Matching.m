addpath("data","dynamics", "plot", "../params", "../../casadi-3.6.7-windows64-matlab2018b")
clear
%%
%Define Dynamic parameters
params = define_params();
m0 = params{1}; I0 = params{2}; mu = params{3}; r= params{4}; d= params{5}; 
thrust_limit= params{6}; kt = params{7}; c_1 = params{8}; c_2 = params{9}; mass_door = params{10}; 
handle_factor = params{11}; inertia = params{12}; r_i_ci = params{13}; n = params{14}; dh = params{15}; gravity = params{16};

m = 5; 
%K = 17; L = 9; core = [9,1];
K = 2 * m - 1; L = m; core = [m, 1];

all_shapes = generate_all_shapes(m, K, L, core);

% NLP parameters
dt = 0.1; N = 100;

x_0 = [0;0;0;0;0;0;0;0];
x_f = [pi/4;0;0;0;0;0;0;0];

% Generate reference trajectory & initial guess of x
t0 = 1; % cosine smoothing
t1 = 1; % hovering time
q_o_ref = generate_q_o_ref(x_0, x_f, N, dt, t0 ,t1);
x_interp = generate_x_interp(x_0, x_f, N, dt, t1);

%%
thrust_scale = 1;
tau_scale = 0;

u_max = thrust_limit * thrust_scale;
u_min = thrust_limit * (- thrust_scale);

tau_min = -0.2 * tau_scale; 
tau_max =  0.2 * tau_scale;

for i= 2:1:min(2,m)
    num_AMs = i;
    shapes = all_shapes{i};
    
    iter = length(shapes);
    all_x_opt = cell(iter, 1);           
    all_u_opt = cell(iter, 1); 
    all_optimal_value = cell(iter, 1);
    all_exit_flag = cell(iter, 1);
    all_processing_time = cell(iter, 1);
    
    nu = 2 + num_AMs*4; zero_us = zeros(nu,1); 
    X_init_guess = [reshape(x_interp',(N+1)*8,1); repmat(zero_us, N, 1)];
    
    % nlp for each shape
    tic;
    for j =1:1:iter
        fprintf('\nnum AMs: %d\nStart shape %d / %d\n', num_AMs, j,iter);
        fprintf('Shape: \n');
        shape = shapes{j};
        disp(shape);
        
        [AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);
        mass =  {mass_door(1), mass_door(2), mass_door(3), AM_mass};
        inertia{4} = AM_inertia;
        r_i_ci{4} = [AM_com(1); 0; AM_com(2)];
    
        %Solve NLP
        [ x_opt, u_opt, optimal_value,exit_flag,processing_time] ....
              = solve_nlp(params, shape, num_AMs, mass, inertia, r_i_ci, q_o_ref, tau_min, tau_max, u_min, u_max, x_0, x_f, X_init_guess, dt, N);
        % Save results for this iteration
        all_x_opt{j} = x_opt;
        all_u_opt{j} = u_opt;
        all_optimal_value{j} = optimal_value;
        all_exit_flag{j} = exit_flag;
        all_processing_time{j} = processing_time;
        fprintf('exit flag: %d\n', exit_flag);
        fprintf('optimal_value: %f\n', optimal_value);
    end

    elapsed_time = toc;  
    fprintf('Total time: %f seconds\n', elapsed_time);
    
    % Save the result
    mkdir data\old_result\ref_1
    file_name = sprintf("data/old_result/ref_1/%d_%d_%d", num_AMs, thrust_scale, tau_scale);
    file_name = sprintf("data/old_result/test.mat");
    save(file_name);
    
    mkdir ..\3D_ver2\data\old_result\ref_1
    file_name = sprintf("../3D_ver2/data/old_result/ref_1/%d_%d_%d", num_AMs,thrust_scale,tau_scale);
    file_name = sprintf("../3D_ver2/data/old_result/test.mat");
    save(file_name);
end
%%
plot_graph

% Video
slow_factor = 1; force_scale = 0.5; do_view = 0 ;
robot = generate_door_ver2(n,dh,r_i_ci, d, gravity, shape, mass,inertia, do_view,q);
%save_plot_tree(robot,dh, params, x_opt,u_opt, dt,N,slow_factor, force_scale, shape)
%plot_tree(robot, dh, params, x_opt,u_opt, dt,N,slow_factor, force_scale, shape)