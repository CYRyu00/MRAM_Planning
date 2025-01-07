num_AMs = 8;
max_serial_num = num_AMs; %12;
max_parallel_num = num_AMs; % 12;
mass_scale = 1;
%file_name = sprintf("Planning_results/%d_%d_%d_optimization_results.mat", num_AMs,max_serial_num,max_parallel_num);
file_name = sprintf("mass_scale_2/7_2_1_optimization_results.mat");

%load(file_name);
%% Assume shapes is a cell array of combinations
% Get the length of each combination
shape_lengths = cellfun(@length, shapes);

% Sort the shapes based on their lengths
[sorted_shape_lengths, sort_indices] = sort(shape_lengths,"descend");  % Get indices of sorted lengths
sort_indices = flip(sort_indices);
sorted_shape_lengths = flip(sorted_shape_lengths);

sorted_shapes = shapes(sort_indices);     % Apply sorting to shapes
sorted_optimal_value = cell2mat(all_optimal_value(sort_indices));
sorted_exit_flag = cell2mat(all_exit_flag(sort_indices));
% Display the sorted shapes

disp('Sorted Shapes by Length:');
disp(sorted_shapes);
%%
close all
figure('Position',[100 100 800 600]);
hold on;
% Separate indices based on the value of sorted_exit_flag
blue_indices = sorted_exit_flag == 1;   % Logical array for blue points
red_indices = sorted_exit_flag ~= 1;    % Logical array for red points

plot(find(blue_indices), sorted_optimal_value(blue_indices), 'bo');
plot(find(red_indices), sorted_optimal_value(red_indices), 'ro');
%plot(sorted_shape_lengths*10,'black')

change_indices = find([true, diff(sorted_shape_lengths) ~= 0]);
% Add vertical lines at the change indices
for idx = change_indices
    xline(idx-0.5, '--k', 'LineWidth', 0.5); % Dashed black vertical line
    if (idx==1)
        txt_disp = 1;
    elseif(idx==2)
        txt_disp = -0.5;
    else 
        txt_disp = 0;
    end
     text(idx-txt_disp, 40, ... %sorted_shape_lengths(idx) * 10, ...
        sprintf('%d', sorted_shape_lengths(idx)), ... % Text is the shape length
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 12, ...
        'Color', 'b');
end

% Add labels and legend
xlabel('Shape');
%xlim([0, length(sorted_exit_flag)+2]);
%ylim([40, 130])
ylabel('Optimal Value');
title('Object function Values for each Shape');
legend('Feasible','Infeasible','Length');
hold off;
%axis tight
%% Find best and worst Shape
% Convert cell arrays to numeric arrays (assuming they are numeric)
all_optimal_value_array = cell2mat(all_optimal_value);
all_exit_flag_array = cell2mat(all_exit_flag);  % Assuming exit_flag is numeric
valid_indices = find(all_exit_flag_array == 1);

if ~isempty(valid_indices)
    valid_optimal_values = all_optimal_value_array(valid_indices);
    [min_value, local_min_index] = min(valid_optimal_values);
    global_min_index = valid_indices(local_min_index);
    
    [max_value, local_max_index] = max(valid_optimal_values);
    global_max_index = valid_indices(local_max_index);

    best_num_up = cell2mat(shapes(global_min_index));
    worst_num_up = cell2mat(shapes(global_max_index));
   
    fprintf('Best Shape (Min) is %.4f at global index %d\n', min_value, global_min_index);
    fprintf('Best Shape: ');
    disp(best_num_up);
    
    fprintf('Worst Shape (Max) is %.4f at global index %d\n', max_value, global_max_index);
    fprintf('Worst Shape: ');
    disp(worst_num_up);
else
    fprintf('No valid exit_flag == 1 found.\n');
end