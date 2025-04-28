%Define Dynamcis and 
for i=1:1:8
object_mass_scale = 2;
AM_mass_scale = 1;
thrust_scale = 1;
define_global_params(object_mass_scale);
[m1, m2, lp, lg, m0, I0,mu,r,d,g,c_cart,c_pole,thrust_limit] = get_global_params();
scale_d = 3;
d = d*scale_d;

num_AMs = i;
max_serial_num = num_AMs; %10;
max_parallel_num = num_AMs;%10;
L_arr = (lp+lg):d:(lp+lg + (max_serial_num -1)*d);
shapes = make_configurations(num_AMs,max_serial_num,max_parallel_num);
%disp(length(shapes)) 2^n-1

x_0 = [0;0;0;0];
x_f = [3;pi/3;0;0];
u_max = thrust_limit *thrust_scale;
u_min = thrust_limit *(-0);
m0=m0*AM_mass_scale;
I0=I0*AM_mass_scale;
dt = 0.05;
N = 100;
%%
iter = length(shapes);
% Preallocate cell arrays or structures to store results
all_x_opt = cell(iter, 1);  % Store x_opt for each iteration
all_u_opt = cell(iter, 1);  % Store u_opt for each iteration
all_optimal_value = cell(iter, 1);  % Store optimal_value for each iteration
all_exit_flag = cell(iter, 1);  % Store exit_flag for each iteration
all_processing_time = cell(iter, 1);  % Store processing_time for each iteration

% setting for initial guess
x_interp = zeros(N+1, 4);
vel_x = (x_f(1) - x_0(1))/(N*dt) ;
vel_theta = (x_f(2) - x_0(2))/(N*dt);
for j = 1:4
    x_interp(:, j) = linspace(x_0(j), x_f(j), N+1)';
end
for j = 1:(N+1)
    x_interp(j, 3:4) = [vel_x;vel_theta];
end
load('1AM_solution.mat')
zero_xs = zeros(length(4),1);

%% nlp for each shape
tic;
for i =1:iter
    num_up = cell2mat(shapes(i));
    fprintf('Start shape %d / %d\n', i,iter);
    fprintf('Shape: ');
    disp(num_up);
    N_u = length(num_up)*2;
    % Make Initial guess of solution
    
    zero_us = zeros(N_u,1);
    
    X_init_guess = [reshape(x_interp',(N+1)*4,1);repmat(zero_us, N, 1)];  % 61.510837 seconds
    
    %u_guess = repmat(u_opt_1AM/num_AMs,1,length(num_up));
    %X_init_guess = [reshape(x_opt_1AM',(N+1)*4,1);reshape(u_guess',N*N_u,1)]; %58.197548  seconds
    
    %X_init_guess = [reshape(x_opt_1AM',(N+1)*4,1);repmat(zero_us, N, 1)]; %57.201678  seconds
    %if i == 1 
    %    X_init_guess = [reshape(x_opt_1AM',(N+1)*4,1);repmat(zero_us, N, 1)];
    %else
    %    X_init_guess = [reshape(x_opt',(N+1)*4,1);repmat(zero_us, N, 1)]; %55.101723  sec
    %end

    %Solve NLP
    [x_opt, u_opt, optimal_value,exit_flag,processing_time] ...
        = solve_nlp(num_up,L_arr,u_min,u_max,x_0,x_f,X_init_guess,dt,N);

    % Save results for this iteration
    all_x_opt{i} = x_opt;
    all_u_opt{i} = u_opt;
    all_optimal_value{i} = optimal_value;
    all_exit_flag{i} = exit_flag;
    all_processing_time{i} = processing_time;
    fprintf('exit flag: %d\n', exit_flag);
end
elapsed_time = toc;  
fprintf('Total time: %f seconds\n', elapsed_time);
%% Save the result
%save('10_10_optimization_results.mat', 'all_x_opt', 'all_u_opt', 'all_optimal_value', 'all_exit_flag', 'all_processing_time');
file_name = sprintf("mass_scale_2/scale_d=3_%d_%d_%d_optimization_results.mat", num_AMs,object_mass_scale, AM_mass_scale);
%elapsed_time = sum(cell2mat(all_processing_time));
global params
save(file_name, 'all_x_opt', 'all_u_opt', 'all_optimal_value', 'all_exit_flag', 'all_processing_time','elapsed_time','shapes'...
    ,'params','L_arr','x_0','x_f','u_min','u_max','dt','N','object_mass_scale','AM_mass_scale','thrust_scale');

end
%% Find best and worst Shape
% Convert cell arrays to numeric arrays (assuming they are numeric)
all_optimal_value_array = cell2mat(all_optimal_value);
all_exit_flag_array = cell2mat(all_exit_flag);  % Assuming exit_flag is numeric

% Find the indices where exit_flag == 1
valid_indices = find(all_exit_flag_array == 1);

% If there are valid indices (exit_flag == 1)
if ~isempty(valid_indices)
    % Filter the optimal values corresponding to exit_flag == 1
    valid_optimal_values = all_optimal_value_array(valid_indices);
    
    % Find the minimum value and its index in the filtered values
    [min_value, local_min_index] = min(valid_optimal_values);
    
    % Convert local_min_index to the global index
    global_min_index = valid_indices(local_min_index);
    
    % Find the maximum value and its index in the filtered values
    [max_value, local_max_index] = max(valid_optimal_values);
    
    % Convert local_max_index to the global index
    global_max_index = valid_indices(local_max_index);
   
    % Extract the corresponding best and worst num_up
    best_num_up = cell2mat(shapes(global_min_index));
    worst_num_up = cell2mat(shapes(global_max_index));
    
    % Display the results
    fprintf('Best Shape (Min) is %.4f at global index %d\n', min_value, global_min_index);
    fprintf('Best Shape: ');
    disp(best_num_up);
    
    fprintf('Worst Shape (Max) is %.4f at global index %d\n', max_value, global_max_index);
    fprintf('Worst Shape: ');
    disp(worst_num_up);
else
    fprintf('No valid exit_flag == 1 found.\n');
end

%% Plot
index = global_min_index;
x_opt = cell2mat(all_x_opt(index));
u_opt = cell2mat(all_u_opt(index));
num_up = cell2mat(shapes(index));

close all
figure 
subplot(2,1,1)
plot(all_optimal_value_array,'.-')
title("object function value")
subplot(2,1,2)
plot(all_exit_flag_array,'*')
title("exit flag")
%% plot
plot_multi_results(x_opt, u_opt,dt,N,L_arr,num_up);
%% Video
slow_factor = 1;
force_scale = 2;

plot_tree_multi(x_opt,u_opt,dt,N,slow_factor,force_scale,L_arr,num_up)