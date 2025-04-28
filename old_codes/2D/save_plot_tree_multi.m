% (x_opt,u_opt, dt,N,slow_factor,scale, L,num_up)
num_AMs = 8;
max_serial_num = num_AMs; % 12;
max_parallel_num = num_AMs; % 12;
file_name = sprintf("Planning_results/%d_%d_%d_optimization_results.mat", num_AMs,max_serial_num,max_parallel_num);
load(file_name);
% Find best and worst Shape
all_optimal_value_array = cell2mat(all_optimal_value);
all_exit_flag_array = cell2mat(all_exit_flag);
valid_indices = find(all_exit_flag_array == 1);
if ~isempty(valid_indices)
    valid_optimal_values = all_optimal_value_array(valid_indices);
    [min_value, local_min_index] = min(valid_optimal_values);
    global_min_index = valid_indices(local_min_index);
    [max_value, local_max_index] = max(valid_optimal_values);

    global_max_index = valid_indices(local_max_index);

    best_num_up = cell2mat(shapes(global_min_index));
    worst_num_up = cell2mat(shapes(global_max_index));
else
    fprintf('No valid exit_flag == 1 found.\n');
end

[m1, m2, lp, lg, m0, I0,mu,r,d,g,c_cart,c_pole,thrust_limit] = get_global_params();

index = global_min_index;
x_opt = cell2mat(all_x_opt(index));
u_opt = cell2mat(all_u_opt(index));
num_up = cell2mat(shapes(index));
% compare best and worst
index = global_min_index;
x_opt1 = cell2mat(all_x_opt(index));
u_opt1 = cell2mat(all_u_opt(index));
num_up1 = cell2mat(shapes(index));
index = global_max_index;
x_opt2 = cell2mat(all_x_opt(index));
u_opt2 = cell2mat(all_u_opt(index));
num_up2 = cell2mat(shapes(index));

taus1 = zeros(length(num_up1),N,2);
for i=1:length(num_up1)
    for j= 1:N
        tau_ = get_tau_distance(x_opt1(j,:)',u_opt1(j,2*i-1:2*i)',L(i));
        taus1(i,j,:)= tau_';
    end
end
sum_taus1 = sum(taus1, 1);

taus2 = zeros(length(num_up2),N,2);
for i=1:length(num_up2)
    for j= 1:N
        tau_ = get_tau_distance(x_opt1(j,:)',u_opt2(j,2*i-1:2*i)',L(i));
        taus2(i,j,:)= tau_';
    end
end
sum_taus2 = sum(taus2, 1); 

%%
close all
slow_factor = 1; scale = 1;
L=L_arr; num_up = worst_num_up;
x = x_opt2;
u = u_opt2;

robot = generate_model_multi_AMs(L,num_up);

figure('Position',[100,100,1000,600])
%subplot(2,2,4)
x0=[0,-pi/2,0,0]';
show(robot,x(1,1:2)'+x0(1:2));

view(2)
ax = gca;
ax.Projection = 'orthographic';
%hold on

global params
mu = params(7);
r = params(8);
d = params(9);
A = [r r -r -r;-r r r -r;mu -mu mu -mu; zeros(2,4);ones(1,4)];

ax.View =[30,20];
axis([-0.5 max(x_opt(:,1))+0.5+d*length(num_up) -1 1 -1 1])


video_filename = 'images/Worst_graph.avi';
v = VideoWriter(video_filename); % Create a video writer object
framesPerSecond = 1/dt*slow_factor;
%rate = rateControl(framesPerSecond);
v.FrameRate = framesPerSecond; % Set frame rate (adjust as needed)
open(v); % Open the video writer
t = (0:1:N-1)*dt;

% Main simulation loop
for i = 1:N
    % Update the robot visualization
    if true 
    subplot(2,2,1)
    hold on
    plot(t, sum_taus1(1,:,1),'b.-');
    plot(t(1:i), sum_taus2(1,1:i,1),'r.-');
    %legend("Best")
    legend("Best","Worst");
    title("Generalized Force")
    ylabel("tau1 [N]")
    xlabel("time [sec]")
    xlim([0 5])
    ylim([-3 4])
    grid on
    hold off
    
    
    subplot(2,2,2)
    hold on
    plot(t, sum_taus1(1,:,2),'b.-');
    plot(t(1:i), sum_taus2(1,1:i,2),'r.-');
    xlim([0 5])
    ylim([-1 3])
    ylabel("tau2 [N·M]")
    xlabel("time [sec]")
    %legend("Best")
    legend("Best","Worst");
    grid on
    hold off
    end
    

    if false
    %subplot(2,2,4)
    show(robot, x(i,1:2)' + x0(1:2), 'PreservePlot', false);
    ax.View = [30,20];%[30,20];
    axis([-0.5 max(x_opt(:,1))+0.5+d*length(num_up) -1 1 -1 1])
    drawnow;
    hold on;

    colors = ['r', 'g', 'b'];
    for j = 1:length(num_up)
        % Get transformation matrix
        tform1 = getTransform(robot, x(i,1:2)' + x0(1:2), sprintf('AM%d-1', j));
        pos1 = tform1(1:3, 4);
        R1 = tform1(1:3, 1:3);

        % Thrust 1
        F_b = A * [u(i, 2*j-1), 0, 0, 0]';
        f_w = R1 * F_b(4:6);
        p1 = pos1 + R1 * [r; r; 0];
        p2 = p1 + f_w * scale;
        plot3([p1(1), p2(1)], [p1(2), p2(2)], [p1(3), p2(3)], colors(mod(j-1, 3) + 1), 'LineWidth', 1);

        % Thrust 2
        F_b = A * [0, u(i, 2*j), 0, 0]';
        f_w = R1 * F_b(4:6);
        p1 = pos1 + R1 * [-r; r; 0];
        p2 = p1 + f_w * scale;
        plot3([p1(1), p2(1)], [p1(2), p2(2)], [p1(3), p2(3)], colors(mod(j-1, 3) + 1), 'LineWidth', 1);
    end
    end

    % Capture the frame
    frame = getframe(gcf);
    writeVideo(v, frame);

    % Pause for the rate if needed
    waitfor(rate);
    hold off
end

% Close video writer
close(v);
disp(['Video saved to ', video_filename]);
hold off;
