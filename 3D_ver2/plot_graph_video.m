addpath("dynamics\","params\" ,"plot\")
%% plot
clear all
load("..\3D_ver2\data\result_9_5\hover\max_iter_1000\10_3_0.mat")
 
close all
figure('Position',[900,100,900,800])
time = (1:1:N+1)*dt;
subplot(4,1,1)
plot(time, x_opt(:,1:4))
hold on
plot(time, qo_desired, "--")
legend({'$q_1$','$q_2$','$q_3$','$q_4$','$q_{1,\mathrm{ref}}$','$q_{2,\mathrm{ref}}$'}, ...
       'Interpreter','latex','FontSize',14);
title("states")
axis tight;

subplot(4,1,2)
plot(time(1:end-1),u_opt(:,1:2))
legend({'$u_1$','$u_2$'}, ...
       'Interpreter','latex','FontSize',14);
title("motor inputs")
axis tight;

[AM_com, AM_mass, AM_inertia]  = get_inertia_double(lau_opt,K,L, core ,m0, I0, d);
mass =  {mass_door(1), mass_door(2), mass_door(1), AM_mass};
inertia = {eye(3)*1, eye(3)*0.1, eye(3)*0.1, AM_inertia, zeros(3,3)};
r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};

wrench = zeros(N,6);
tau = zeros(N,n);
for i=1:N
    wrench(i,:) = map_u2wrench_double( u_opt(i,3:end)',lau_opt,K,L, core , mu , r , d);
    q = x_opt(i,1:4)'; qd = x_opt(i,5:8)'; qdd = (x_opt(i+1,5:8) -x_opt(i+1,5:8) )'/dt; 
    F_ext = wrench(i,:)';
    tau(i,:) =  [-c_1*qd(1);(-c_2*qd(2) -kt*q(2) + mass{2}*handle_factor); u_opt(i,1); u_opt(i,2)] + ...
        newton_euler_inverse_dynamics_double(n, dh, mass, inertia, r_i_ci, gravity, q, qd, qdd, F_ext);
end
subplot(4,1,3)
plot(time(1:end-1),wrench)
legend({'$m_x$','$m_y$','$m_z$','$f_x$','$f_y$','$f_z$'}, ...
       'Interpreter','latex','FontSize',14);

title("Wrench wrt. AM frame")
axis tight;

subplot(4,1,4)
plot(time(1:end-1),tau)
axis tight
legend({'$\tau_1$','$\tau_2$','$\tau_3$','$\tau_4$'}, ...
       'Interpreter','latex','FontSize',14);
title("generalized force")
axis tight;

%plot 3d video
do_view=1; q =  [0;0;0;0]; g=[0;0;-9.81];
robot = generate_door(n,dh,r_i_ci,d, g, lau_opt, core, mass,inertia, do_view,q);

slow_factor =1; force_scale = 0.2;
save_plot_tree(robot,dh, params, x_opt,u_opt, dt,N,slow_factor, force_scale, lau_opt, core, K, L)
%plot_tree(robot,dh, params, x_opt,u_opt, dt,N,slow_factor, force_scale, lau_opt, core, K, L)
