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

% Plot f* wrt inertia
close all
figure 
hold on
for i = 1:length(all_exit_flag_array)
    if all_exit_flag_array(i) == 1
        plot(i, all_optimal_value_array(i), 'b.') % Red for exit_flag == 1
    else
        plot(i, all_optimal_value_array(i), 'r.') % Blue otherwise
    end
end
hold off

Ixx_array = [];
Iyy_array = [];
Izz_array = [];
for i = 1:length( all_optimal_value_array)
    [tmp, tmp2, Inertia] = get_inertia_ver2(shapes{i} ,m0, I0, d);
    Ixx_array = [Ixx_array;Inertia(1,1)];
    Iyy_array = [Iyy_array;Inertia(2,2)];
    Izz_array = [Izz_array;Inertia(3,3)];
end

%ascending order
[sorted_Ixx, sortIdx_xx] = sort(Ixx_array);
sorted_xx_all_optimal_value_array = all_optimal_value_array(sortIdx_xx);
[sorted_Iyy, sortIdx_yy] = sort(Iyy_array);
sorted_yy_all_optimal_value_array = all_optimal_value_array(sortIdx_yy);
[sorted_Izz, sortIdx_zz] = sort(Izz_array);
sorted_zz_all_optimal_value_array = all_optimal_value_array(sortIdx_zz);

figure 
subplot(3,1,1)
hold on
for i = 1:length(all_exit_flag_array)
    if all_exit_flag_array(i) == 1
        plot(sorted_Ixx(i),sorted_xx_all_optimal_value_array(i), 'b.') % Red for exit_flag == 1
    else
        plot(sorted_Ixx(i),sorted_xx_all_optimal_value_array(i), 'r.') % Blue otherwise
    end
end
hold off
title("I_{xx} - Cost func. value")

subplot(3,1,2)
hold on
for i = 1:length(all_exit_flag_array)
    if all_exit_flag_array(i) == 1
        plot(sorted_Iyy(i), sorted_yy_all_optimal_value_array(i), 'b.') % Red for exit_flag == 1
    else
        plot(sorted_Iyy(i), sorted_yy_all_optimal_value_array(i), 'r.') % Blue otherwise
    end
end
hold off
title("I_{yy} - Cost func. value")

subplot(3,1,3)
hold on
for i = 1:length(all_exit_flag_array)
    if all_exit_flag_array(i) == 1
        plot(sorted_Izz(i), sorted_zz_all_optimal_value_array(i), 'b.') % Red for exit_flag == 1
    else
        plot(sorted_Izz(i), sorted_zz_all_optimal_value_array(i), 'r.') % Blue otherwise
    end
end
hold off
title("I_{zz} - Cost func. value")
%% plot
index = global_min_index;
%index = 369;
disp(shapes{index})
x_opt = cell2mat(all_x_opt(index));
u_opt = cell2mat(all_u_opt(index));
shape = cell2mat(shapes(index));

[AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);
mass =  {mass_door(1), mass_door(2), mass_door(3), AM_mass};
inertia = {eye(3)*1, eye(3)*0.1, eye(3)*0.1, AM_inertia, zeros(3,3)};
r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};

figure('Position',[800,100,800,600])
time = (1:1:N+1)*dt;
subplot(4,1,1)
plot(time, x_opt(:,1:4))
hold on
plot(time, q_o_ref, "--")
legend({'$q_1$','$q_2$','$q_3$','$q_4$','$q_{1,\mathrm{ref}}$','$q_{2,\mathrm{ref}}$','$q_{3,\mathrm{ref}}$', '$q_{4,\mathrm{ref}}$'}, ...
       'Interpreter','latex','FontSize',12);
title("states")
axis tight;

subplot(4,1,2)
plot(time(1:end-1),u_opt(:,1:2))
legend({'$u_1$','$u_2$'}, ...
       'Interpreter','latex','FontSize',12);
title("motor inputs")
axis tight;

wrench = zeros(N,6);
tau = zeros(N,n);
for i=1:N
    wrench(i,:) = map_u2wrench_double( u_opt(i,3:end)', shape , mu , r , d);
    q = x_opt(i,1:4)'; qd = x_opt(i,5:8)'; qdd = (x_opt(i+1,5:8) -x_opt(i+1,5:8) )'/dt; 
    F_ext = wrench(i,:)';
    
    tau(i,:) =  [0; -kt*q(2) + mass{2}*handle_factor; u_opt(i,1); u_opt(i,2)] + ...
        newton_euler_inverse_dynamics_double(n, dh, mass, inertia, r_i_ci, gravity, q, qd, qdd, F_ext);
end
subplot(4,1,3)
plot(time(1:end-1),wrench)
legend({'$m_x$','$m_y$','$m_z$','$f_x$','$f_y$','$f_z$'}, ...
       'Interpreter','latex','FontSize',12);
title("Wrench exp. AM frame")
axis tight;

subplot(4,1,4)
plot(time(1:end-1),tau)
axis tight
legend({'$\tau_1$','$\tau_2$','$\tau_3$','$\tau_4$'}, ...
       'Interpreter','latex','FontSize',12);
title("generalized force")
axis tight;


