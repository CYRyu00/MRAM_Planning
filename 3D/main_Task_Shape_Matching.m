%Define Dynamic parameters and shapes
n = 4;
% dynamic parameters of each module 
params = define_params();
m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5}; thrust_limit= params{6};kt=params{7};c_1=params{8};c_2=params{9};

dh = [0,0,0.95,0;   % [alpha, a, d, theta]
      -pi/2, 0.9 , 0,0;
      0,-0.1,0.23,pi/2;
      pi/2,0,0,-pi/2;
      pi/2,0,0,0];
gravity = [0;0;-9.81];

m=6;
all_shapes = generate_all_shapes(m);

% Show example rigidbodytree
shape = all_shapes{6}{10};
[AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);
mass =  {10, 1, 0.5, AM_mass};
inertia = {eye(3)*1, eye(3)*0.01, eye(3)*0.01, AM_inertia, zeros(3,3)};
r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};

do_view =0; q=[pi/6;-pi/3;pi/3;0];
robot = generate_door_ver2(n,dh,r_i_ci, d, gravity, shape, mass,inertia, do_view,q);

% NLP parameters
dt = 0.1;
N = 100;

x_0 = [0;0;0;0;0;0;0;0];
x_f = [pi/4;0;0;0;0;0;0;0];
qo_desired = zeros(2,N+1);
T = N*dt; t0 = 1; %1sec
for j=1:N+1
    t = (j-1)*dt;
    if t < t0 
        qo_desired(:,j) = [0 ; pi/6*(cos(pi/t0*t) - 1)/2];
    elseif t < T-t0
        qo_desired(:,j) = [(x_f(1)-x_0(1))/(T-t0)*(t-t0) ;-pi/6];
    else
        qo_desired(:,j) = [(x_f(1)-x_0(1))/(T-t0)*(t-t0) ;pi/6*(cos(pi/t0*t -pi/t0*T) - 1)/2];
    end
end

thrust_scale = 1;
u_max = thrust_limit *thrust_scale;
u_min = thrust_limit *(-thrust_scale);
tau_scale = 1;
tau_min = -0.2 *tau_scale ; 
tau_max =  0.2 *tau_scale ;

%%
for i=6:1:6
num_AMs = i;
shapes = all_shapes{i};

%
iter = length(shapes);
all_x_opt = cell(iter, 1);  % Store x_opt for each iteration
all_u_opt = cell(iter, 1);  % Store u_opt for each iteration
all_optimal_value = cell(iter, 1);  % Store optimal_value for each iteration
all_exit_flag = cell(iter, 1);  % Store exit_flag for each iteration
all_processing_time = cell(iter, 1);  % Store processing_time for each iteration

x_interp = zeros(N+1, 8);
vel_x1 = (x_f(1) - x_0(1))/(N*dt); vel_x2 = (x_f(2) - x_0(2))/(N*dt);
vel_x3 = (x_f(3) - x_0(3))/(N*dt); vel_x4 = (x_f(4) - x_0(4))/(N*dt);
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
for j =1:1:iter
    fprintf('Start shape %d / %d\n', j,iter);
    fprintf('Shape: \n');
    shape = shapes{j};
    disp(shape);
    
    [AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);
    mass =  {10, 1, 0.5, AM_mass};
    inertia = {eye(3)*1, eye(3)*0.1, eye(3)*0.1, AM_inertia, zeros(3,3)};
    r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};

    %Solve NLP
    [ x_opt, u_opt, optimal_value,exit_flag,processing_time] ....
          = solve_nlp(params,shape,num_AMs,dh, gravity,mass,inertia,r_i_ci, ...
                      qo_desired, tau_min, tau_max,u_min,u_max,x_0,x_f,X_init_guess,dt,N);
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
file_name = sprintf("result/%d_%d_%d_optimization_results.mat", num_AMs,thrust_scale,tau_scale);
elapsed_time = sum(cell2mat(all_processing_time));

save(file_name, 'all_x_opt', 'all_u_opt', 'all_optimal_value', 'all_exit_flag', 'all_processing_time','elapsed_time','shapes'...
    ,'params','x_0','x_f','u_min','u_max','tau_min', 'tau_max','dt','N','thrust_scale');

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
    fprintf('Best Shape: \n');
    disp(best_shape);
    
    fprintf('Worst Shape (Max) is %.4f at global index %d\n', max_value, global_max_index);
    fprintf('Worst Shape: \n');
    disp(worst_shape);
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
index = global_max_index;
index = 369;
disp(shapes{index})
x_opt = cell2mat(all_x_opt(index));
u_opt = cell2mat(all_u_opt(index));
shape = cell2mat(shapes(index));

[AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);
mass =  {10, 1, 0.5, AM_mass};
inertia = {eye(3)*1, eye(3)*0.1, eye(3)*0.1, AM_inertia, zeros(3,3)};
r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};

figure('Position',[800,100,800,600])
time = (1:1:N+1)*dt;
subplot(4,1,1)
plot(time, x_opt(:,1:4))
hold on
plot(time, qo_desired, "--")
legend("q_1","q_2","q_3","q_4","q_{1,desired}","q_{2,desired}")
title("states")
axis tight;

subplot(4,1,2)
plot(time(1:end-1),u_opt(:,1:2))
legend("u_1","u_2")
title("motor inputs")
axis tight;

wrench = zeros(N,6);
tau = zeros(N,n);
for i=1:N
    wrench(i,:) = map_u2wrench_double( u_opt(i,3:end)', shape , mu , r , d);
    q = x_opt(i,1:4)'; qd = x_opt(i,5:8)'; qdd = (x_opt(i+1,5:8) -x_opt(i+1,5:8) )'/dt; 
    F_ext = wrench(i,:)';
    
    tau(i,:) =  [0; -kt*q(2) + mass{2}*0.80*9.81*0.040; u_opt(i,1); u_opt(i,2)] + ...
        newton_euler_inverse_dynamics_double(n, dh, mass, inertia, r_i_ci, gravity, q, qd, qdd, F_ext);
end
subplot(4,1,3)
plot(time(1:end-1),wrench)
legend("m_x","m_y","m_z","f_x","f_y","f_z")
title("Wrench exp. AM frame")
axis tight;

subplot(4,1,4)
plot(time(1:end-1),tau)
axis tight
legend("tau1","tau2","tau3","tau4")
title("generalized force")
axis tight;


% Video
slow_factor = 1;
force_scale = 0.5;


do_view = 0 ;
robot = generate_door_ver2(n,dh,r_i_ci, d, gravity, shape, mass,inertia, do_view,q);

%save_plot_tree(robot,dh, params, x_opt,u_opt, dt,N,slow_factor, force_scale, shape)
plot_tree(robot, dh, params, x_opt,u_opt, dt,N,slow_factor, force_scale, shape)