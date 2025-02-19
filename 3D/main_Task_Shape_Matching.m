%Define Dynamic parameters and shapes
n = 4;
% dynamic parameters of each module 
params = define_params();
m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5};
thrust_limit= params{6};

n=4;
dh = [0,0,0.95,0;   % [alpha, a, d, theta]
      -pi/2, 0.9 , 0,0;
      0,-0.1,0.23,pi/2;
      pi/2,0,0,-pi/2;
      pi/2,0,0,0];
gravity = [0;0;-9.81];

params = define_params();
m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5};thrust_limit= params{6};

m=7;
all_shapes = generate_all_shapes(m);
shape = all_shapes{7}{10};
[AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);
mass =  {10, 1, 0.5, AM_mass};
inertia = {eye(3)*1, eye(3)*0.1, eye(3)*0.1, AM_inertia, zeros(3,3)};
r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};

do_view =1; q=[pi/6;-pi/3;pi/3;-pi/4];
robot = generate_door_ver2(n,dh,r_i_ci, gravity, mass,inertia, do_view,q);

x_0 = [0;0;0;0;0;0;0;0];
x_f = [pi/3;-pi/6;pi/6;0;0;0;0;0];
thrust_scale = 5;
u_max = thrust_limit *thrust_scale;
u_min = thrust_limit *(-thrust_scale);
tau_min = -0.75; tau_max = 0.75;
dt = 0.1;
N = 100;
%%
for i=5:1:5

%thrust_scale = inf;
num_AMs = i;
shapes = all_shapes{i};

%x_0 = [0;0;0;0;0;0;0;0];
%x_f = [pi/3;0;0;0;0;0;0;0];
%u_max = thrust_limit *thrust_scale;
%u_min = thrust_limit *(-thrust_scale);

%dt = 0.05;
%N = 100;
%
iter = length(shapes);
all_x_opt = cell(iter, 1);  % Store x_opt for each iteration
all_u_opt = cell(iter, 1);  % Store u_opt for each iteration
all_optimal_value = cell(iter, 1);  % Store optimal_value for each iteration
all_exit_flag = cell(iter, 1);  % Store exit_flag for each iteration
all_processing_time = cell(iter, 1);  % Store processing_time for each iteration

x_interp = zeros(N+1, 8);
vel_x1 = (x_f(1) - x_0(1))/(N*dt); vel_x2 = (x_f(2) - x_0(2))/(N*dt);
vel_x3 = (x_f(3) - x_0(3))/(N*dt);vel_x4 = (x_f(4) - x_0(4))/(N*dt);
for k = 1:8
    x_interp(:, k) = linspace(x_0(k), x_f(k), N+1)';
end
for k = 1:(N+1)
    x_interp(k, 5:8) = [vel_x1;vel_x2;vel_x3;vel_x4];
end

nu = 2 + num_AMs*4; zero_us = zeros(nu,1); 
X_init_guess = [reshape(x_interp',(N+1)*8,1);repmat(zero_us, N, 1)];

%% nlp for each shape
tic;
for j =1:iter
    fprintf('Start shape %d / %d\n', j,iter);
    fprintf('Shape: \n');
    shape = shapes{j};
    disp(shape);
    
    mass =  {10, 1, 0.5, AM_mass};
    inertia = {eye(3)*1, eye(3)*0.1, eye(3)*0.1, AM_inertia, zeros(3,3)};
    r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};

    
    %u_guess = repmat(u_opt_1AM/num_AMs,1,length(num_up));
    %X_init_guess = [reshape(x_opt_1AM',(N+1)*4,1);reshape(u_guess',N*N_u,1)]; %58.197548  seconds
    %X_init_guess = [reshape(x_opt_1AM',(N+1)*4,1);repmat(zero_us, N, 1)]; %57.201678  seconds
    %if i == 1 
    %    X_init_guess = [reshape(x_opt_1AM',(N+1)*4,1);repmat(zero_us, N, 1)];
    %else
    %    X_init_guess = [reshape(x_opt',(N+1)*4,1);repmat(zero_us, N, 1)]; %55.101723  sec
    %end

    %Solve NLP
    [x_opt, u_opt, optimal_value,exit_flag,processing_time] ....
          = solve_nlp(params,shape,num_AMs,dh, gravity,mass,inertia,r_i_ci, ...
                      tau_min, tau_max,u_min,u_max,x_0,x_f,X_init_guess,dt,N);
    % Save results for this iteration
    all_x_opt{j} = x_opt;
    all_u_opt{j} = u_opt;
    all_optimal_value{j} = optimal_value;
    all_exit_flag{j} = exit_flag;
    all_processing_time{j} = processing_time;
    fprintf('exit flag: %d\n', exit_flag);
end
elapsed_time = toc;  
fprintf('Total time: %f seconds\n', elapsed_time);
%% Save the result
%save('10_10_optimization_results.mat', 'all_x_opt', 'all_u_opt', 'all_optimal_value', 'all_exit_flag', 'all_processing_time');
%file_name = sprintf("mass_scale_2/scale_d=3_%d_%d_%d_optimization_results.mat", num_AMs,object_mass_scale, AM_mass_scale);
%elapsed_time = sum(cell2mat(all_processing_time));
%global params
%save(file_name, 'all_x_opt', 'all_u_opt', 'all_optimal_value', 'all_exit_flag', 'all_processing_time','elapsed_time','shapes'...
%    ,'params','L_arr','x_0','x_f','u_min','u_max','dt','N','object_mass_scale','AM_mass_scale','thrust_scale');

end
%% Find best and worst Shape
all_optimal_value_array = cell2mat(all_optimal_value);
all_exit_flag_array = cell2mat(all_exit_flag); 

valid_indices = find(all_exit_flag_array == 1);
if ~isempty(valid_indices)

    valid_optimal_values = all_optimal_value_array(valid_indices);
    [min_value, local_min_index] = min(valid_optimal_values);
    global_min_index = valid_indices(local_min_index);
    
    [max_value, local_max_index] = max(valid_optimal_values);
    global_max_index = valid_indices(local_max_index);
   
    best_shape = cell2mat(shapes(global_min_index));
    worst_shape = cell2mat(shapes(global_max_index));
    
    % Display the results
    fprintf('Best Shape (Min) is %.4f at global index %d\n', min_value, global_min_index);
    fprintf('Best Shape: ');
    disp(best_num_up);
    
    fprintf('Worst Shape (Max) is %.4f at global index %d\n', max_value, global_max_index);
    fprintf('Worst Shape: ');
    disp(worst_shape);
else
    fprintf('No valid exit_flag == 1 found.\n');
end

%% Plot
index = global_min_index;
x_opt = cell2mat(all_x_opt(index));
u_opt = cell2mat(all_u_opt(index));
shape = cell2mat(shapes(index));

close all
figure 
subplot(2,1,1)
plot(all_optimal_value_array,'.-')
title("object function value")
subplot(2,1,2)
plot(all_exit_flag_array,'*')
title("exit flag")
%% plot
%plot_multi_results(x_opt, u_opt,dt,N,L_arr,num_up);
%% Video
slow_factor = 1;
force_scale = 2;

close all
plot_tree(robot, x_opt,u_opt, dt,N,slow_factor,force_scale, shape)