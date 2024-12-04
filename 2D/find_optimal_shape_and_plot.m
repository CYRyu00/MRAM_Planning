num_AMs = 8;
max_serial_num = num_AMs; % 12;
max_parallel_num = num_AMs; % 12;
file_name = sprintf("Planning_results/%d_%d_%d_optimization_results.mat", num_AMs,max_serial_num,max_parallel_num);
%file_name = sprintf("10_10_10_optimization_results.mat", num_AMs,max_serial_num,max_parallel_num);

load(file_name);

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

close all
figure 
subplot(2,1,1)
plot(all_optimal_value_array,'.-')
title("object function value")
subplot(2,1,2)
plot(all_exit_flag_array,'*')
title("exit flag")
%% plot

[m1, m2, lp, lg, m0, I0,mu,r,d,g,c_cart,c_pole,thrust_limit] = get_global_params();

index = global_min_index;
x_opt = cell2mat(all_x_opt(index));
u_opt = cell2mat(all_u_opt(index));
num_up = cell2mat(shapes(index));

%plot a result4
% plot_multi_results(x_opt, u_opt,dt,N,L_arr,num_up);

% compare best and worst
index = global_min_index;
x_opt1 = cell2mat(all_x_opt(index));
u_opt1 = cell2mat(all_u_opt(index));
num_up1 = cell2mat(shapes(index));

index = global_max_index;
x_opt2 = cell2mat(all_x_opt(index));
u_opt2 = cell2mat(all_u_opt(index));
num_up2 = cell2mat(shapes(index));
plot_compare_results(x_opt1, u_opt1,dt,N,L_arr,num_up1,x_opt2, u_opt2,num_up2);
%% Video
slow_factor = 1;
force_scale = 1;

plot_tree_multi(x_opt,u_opt,dt,N,slow_factor,force_scale,L_arr,num_up)