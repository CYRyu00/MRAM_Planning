%% Mass-scale = 5
num_AMs = [7 8 9 10 11 12];
optimal_val = [61.1877, 50.5693 ,47.3755, 45.3840 ,44.3912, 43.7472 ];
best_shape = {[ 1     1     1     1     1     1     1],
    [ 1     1     1     1     1     1     1     1],
    [ 1     1     1     1     1     1     1     1     1],
    [ 2     1     1     1     1     1     1     1     1],
    [ 2     1     1     1     1     1     1     1     1     1],
    [ 3     1     1     1     1     1     1     1     1     1]
 };

figure(1)
plot(num_AMs,optimal_val,'.-')
%% Mass-scale = 1
num_AMs = 2:1:10;
optimal_val = [33.9781, 21.9980 ,21.9122, 22.6126 ,23.5500, 24.4539, 25.3371 ,26.2600,27.1391 ];
best_shape = {[  1     1],
    [ 1     1     1],
    [1     1     1     1],
    [1     1     1     1     1],
    [2     1     1     1     1],
    [2     1     1     1     1     1],
    [3     1     1     1     1     1],
    [4     1     1     1     1     1],
    [4     1     1     1     1     1     1]
 };

figure(1)
plot(num_AMs,optimal_val,'o--')
%% Mass-scale = 2
num_AMs = 4:1:10;
optimal_val = [30.1502 , 27.7084 ,27.3313 , 27.6444  ,28.1781, 28.8074, 29.4050];
best_shape = {  [1     1     1     1],
                [1     1     1     1     1],
                [1     1     1     1     1     1],
                [1     1     1     1     1     1     1],
                [2     1     1     1     1     1     1],
                [2     1     1     1     1     1     1     1],
                [ 3     1     1     1     1     1     1     1],
 };

figure(1)
plot(num_AMs, optimal_val, 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.5, 0.5, 0.5], ...
    'LineStyle', '--', 'Color', 'b', 'LineWidth', 0.70)

xlim([3 11])
ylim([25 33])
legend("Optimal Value for Task")
xlabel("m: number of AMs")
ylabel("Object Function Value")
title("Object Function Value of Best Shape")

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
