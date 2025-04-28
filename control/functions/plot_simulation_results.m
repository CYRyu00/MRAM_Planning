function plot_simulation_results(t_sim, x_sim, x_d_interp, u_sim, u_d_interp, disturb_sim, shape, dh, gravity, n, nx, nu, N_sim, dt_sim, delta_inertia, delta_k, sigma, mean, max_val)
% plot trajectory and simulated x and u
figure('Position',[1100 520 800 480])
colors = lines(nx);
legend_entries = cell(1, nx/2);
subplot(3,1,1)
hold on
for i = 1:nx/2
    % Plot simulation
    plot(t_sim, x_sim(:,i), 'Color', colors(i,:), 'LineWidth', 1.0);
    legend_entries{2*i-1} = sprintf('$q_%d$', i);
    
    % Plot reference
    plot(t_sim, x_d_interp(:,i), '--', 'Color', colors(i,:), 'LineWidth', 1.0);
    legend_entries{2*i} = sprintf('$q_{%d,d}$', i);
end
legend(legend_entries, 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 10)
xlabel('Time [s]'); ylabel('State Value')
title_name = sprintf(['x_{sim} vs x_{desired}, \nMass X %.2f, Spring x %.2f,\n Disturbance : ' ...
                        'Sigma: %.2f, Mean: %.2f, Max val :%.2f'],delta_inertia, delta_k , sigma, mean ,max_val);
title(title_name , 'FontSize',10)
grid on

subplot(3,1,2)
legend_entries = cell(1, nx/2);
hold on
for i = (nx/2+1):nx
    % Plot simulation
    plot(t_sim, x_sim(:,i), 'Color', colors(i,:), 'LineWidth', 1.0);
    legend_entries{2*i-1 -nx} = sprintf('$\\dot{q}_%d$', i);
    
    % Plot reference
    plot(t_sim, x_d_interp(:,i), '--', 'Color', colors(i,:), 'LineWidth', 1.0);
    legend_entries{2*i -nx} = sprintf('$\\dot{q}_{%d,d}$', i);
end
legend(legend_entries, 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 10)
xlabel('Time [s]'); ylabel('State Value')
grid on

subplot(3,1,3)
colors = lines(nu);
hold on
for i = 1:nu
    plot(t_sim(1:end-1), u_sim(:,i), 'Color', colors(i,:), 'LineWidth', 1.0);    
    plot(t_sim(1:end-1), u_d_interp(:,i), '--', 'Color', colors(i,:), 'LineWidth', 1.0);
end
title("u_{sim} vs u_{desired}")
xlabel('Time [s]'); ylabel('Input Value')
grid on

% plot wrench and generalized force
wrench = zeros(N_sim,6);
tau = zeros(N_sim,n);
params = define_params();
m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5}; thrust_limit= params{6}; kt=params{7};c_1=params{8};c_2=params{9}; mass_door = params{10}; handle_factor = params{11};
kt_double = kt * delta_k;
mass_door = mass_door * delta_inertia;
[AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d); 

mass_double =  {mass_door(1), mass_door(2), mass_door(3), AM_mass};
inertia_double = {eye(3)*1, eye(3)*0.1, eye(3)*0.1, AM_inertia, zeros(3,3)};
r_i_ci_double = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};

for i=1:N_sim
    wrench(i,:) = map_u2wrench_double( u_sim(i,:)', shape , mu , r , d);
    q = x_sim(i,1:4)'; qd = x_sim(i,5:8)'; qdd = (x_sim(i+1,5:8) -x_sim(i+1,5:8) )'/dt_sim; 
    F_ext = wrench(i,:)';
    tau(i,:) =  [-c_1*qd(1); full(-c_2*qd(2) -kt_double*q(2) + mass_double{2}*handle_factor); u_sim(i,1); u_sim(i,2)] + ...
        newton_euler_inverse_dynamics_double(n, dh, mass_double, inertia_double, r_i_ci_double, gravity, q, qd, qdd, F_ext) ...
            + disturb_sim(i,:)';
end
figure('Position',[1100 50 800 400])
subplot(2,1,1)
plot(t_sim(1:end-1),wrench)

legend({'$m_x$', '$m_y$', '$m_z$', '$f_x$', '$f_y$', '$f_z$'}, ...
       'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 10);
title("Wrench exp. AM frame")
axis tight;
grid on

subplot(2,1,2)
colors = lines(n);
hold on
for i = 1:n
    plot(t_sim(1:end-1), tau(:,i), 'Color', colors(i,:), 'LineWidth', 1.0);    
    plot(t_sim(1:end-1), disturb_sim(:,i), '.', 'Color', colors(i,:), 'LineWidth', 1.0);
end
axis tight
legend({'$\tau_1$','$d_1$', '$\tau_2$','$d_2$', '$\tau_3$', '$d_3$','$\tau_4$','$d_4$',}, ...
       'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 10);

title("generalized force")
axis tight;
grid on

end
