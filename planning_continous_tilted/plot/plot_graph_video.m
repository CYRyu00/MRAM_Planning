addpath("dynamics","../params" ,"plot")
% plot
%clear all
%load("..\3D_ver2\data\result_9_5\hover\max_iter_1000\5_3_0.mat")
 
close all
figure('Position',[900,100,900,800])
time = (1:1:N+1)*dt;
subplot(4,1,1)
colors = lines(4);  % You can also use your own RGB values

hold on
for i = 1:4
    plot(time, x_opt(:, i), 'Color', colors(i,:), 'LineWidth', 1.5)        % q_i
    plot(time, q_o_ref(i, :), '--', 'Color', colors(i,:), 'LineWidth', 1)  % q_i,ref
end

legend({'$q_1$','$q_{1,\mathrm{ref}}$', ...
        '$q_2$','$q_{2,\mathrm{ref}}$', ...
        '$q_3$','$q_{3,\mathrm{ref}}$', ...
        '$q_4$','$q_{4,\mathrm{ref}}$'}, ...
       'Interpreter','latex','FontSize',12);

title("states")
axis tight

subplot(4,1,2)
plot(time(1:end-1), u_opt,'LineWidth', 1.5)
title("thrusts")
axis tight;

[AM_com, AM_mass, AM_inertia]  = get_inertia_duo_double(rho_opt,K,L, core ,m0, I0, d);
mass = {mass_door(1), mass_door(2), mass_door(3), AM_mass};
inertia{4} = AM_inertia;
r_i_ci{4} = [AM_com(1); r_i_ci{4}(2); AM_com(2)];

wrench = zeros(N,6);
tau = zeros(N,n);
for i=1:N
    wrench(i,:) = map_u2wrench_duo_double(u_opt(i,:)', rho_opt, K, L, core, mu, r, d, theta);
    q = x_opt(i,1:4)'; qd = x_opt(i,5:8)'; qdd = (x_opt(i+1,5:8) -x_opt(i+1,5:8) )'/dt; 
    F_ext = wrench(i,:)';
    tau(i,:) =  [-c_1*qd(1);(-c_2*qd(2) -kt*q(2) + mass{2}*handle_factor); 0; 0] + ...
        newton_euler_inverse_dynamics_double(n, dh, mass, inertia, r_i_ci, gravity, q, qd, qdd, F_ext);
end
subplot(4,1,3)
plot(time(1:end-1), wrench, 'LineWidth', 1.5)
legend({'$m_x$','$m_y$','$m_z$','$f_x$','$f_y$','$f_z$'}, ...
       'Interpreter','latex','FontSize',12);
title("Wrench wrt. AM frame")
axis tight;

subplot(4,1,4)
plot(time(1:end-1), tau, 'LineWidth', 1.5)
axis tight
legend({'$\tau_1$','$\tau_2$','$\tau_3$','$\tau_4$'}, ...
       'Interpreter','latex','FontSize',12);
title("generalized force")
axis tight;

%plot 3d video
do_view = 1; q = [0;0;0;0]; g = [0; 0; -9.81];
robot = generate_door(n, dh, r_i_ci, d, g, rho_opt, core, mass, inertia, do_view, q);

slow_factor = 1; force_scale = 0.2;
%save_plot_tree(robot, dh, params, theta, x_opt, u_opt, dt, N, slow_factor, force_scale, rho_opt, core, K, L)
plot_tree(robot, dh, params, theta, x_opt, u_opt, dt, N, slow_factor, force_scale, rho_opt, core, K, L)
