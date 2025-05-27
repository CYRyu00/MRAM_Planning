addpath("dynamics", "functions", "../params" )
clear; close all
params = define_params();
m0 = params{1}; I0 = params{2}; mu = params{3}; r= params{4}; d= params{5};
thrust_limit= params{6}; kt = params{7}; c_1 = params{8}; c_2 = params{9}; mass_door = params{10};
handle_factor = params{11}; inertia = params{12}; r_i_ci = params{13}; n = params{14}; dh = params{15}; gravity = params{16};
%%
alpha = 5e0; beta = 1e-3;

w_n = 2 * pi * 0.5; damp = 0.90;
b = 2* m0 * damp * w_n; k = m0 * w_n^2; 

w_n = 2 * pi * 5; damp = 0.1;
b_w = 2* m0 * damp * w_n; k_w = m0 * w_n^2; 

Y_hover = [1; 2; 3] / 10; w3_hover = 5 / 180 *pi; 

L_f = eye(3) * 10;
L_tau = eye(3) * 100;

D = [0.01; 0; 0.03];
eta = 1 / 180 * pi;
%%
N = 100; dt = 0.1;
dt_sim = 0.01;
N_sim = N * dt / dt_sim;
times = [];

s = sin(eta); c = cos(eta);
A_eta = [mu*s/sqrt(2) - r*c, mu*s/sqrt(2) - r*c, - mu*s/sqrt(2) + r*c, - mu*s/sqrt(2) + r*c;
         - mu*s/sqrt(2) + r*c, mu*s/sqrt(2) - r*c, mu*s/sqrt(2) - r*c, - mu*s/sqrt(2) + r*c;
         mu*c + sqrt(2)*r*s, - mu*c - sqrt(2)*r*s, mu*c + sqrt(2)*r*s, - mu*c - sqrt(2)*r*s;
         s/sqrt(2), - s/sqrt(2), - s/sqrt(2), s/sqrt(2);
         - s/sqrt(2), - s/sqrt(2), s/sqrt(2), s/sqrt(2);
         c, c, c, c];
A = [r r -r -r;  -r r r -r;mu -mu mu -mu; 1 1 1 1];


Y = [0 ; 0; 0]; Yd = [0 ; 0; 0];
w = [0 ; 0; 0]; wd = [0 ; 0; 0];
R = eye(3, 3);
f_hat = zeros(3, 1); tau_hat = zeros(3, 1); int_wd_des = zeros(3, 1);

Y_hist = []; Yd_hist = []; w_hist = []; wd_hist = []; R_hist = cell(N_sim, 1); wd_des_hist = [];
f_hat_hist = []; tau_hat_hist = []; f_e_hist = []; tau_e_hist = [];
thrusts_hist = [];

[Y_des, Yd_des, Ydd_des, w3_des, w3d_des] = get_traj(Y_hover, w3_hover, N_sim, dt_sim);

for i = 1:N_sim/1
    % External wrench & wrench estimator
    % TODO estimation by observed w_dot, w or noisy
    f_e = [0; 0; 0.05] * sin(i/100); tau_e = [3; 2; 1] * 1e-4 * sin(i/100);
    f_hat  = f_hat  + L_f * (f_e  - f_hat) * dt_sim;
    tau_hat = tau_hat + L_tau * (tau_e - tau_hat) * dt_sim;
    
    f_hat_hist = [f_hat_hist, f_hat];
    tau_hat_hist = [tau_hat_hist, tau_hat];
    f_e_hist = [f_e_hist, f_e];
    tau_e_hist = [tau_e_hist, tau_e];
    
    % Tool position and w3 control
    u_Y = m0 * Ydd_des(:, i) - b * (Yd - Yd_des(:, i)) - k * (Y - Y_des(:, i)) - f_hat;
    nu_ = R' * (u_Y + m0 * norm(gravity) * [0; 0; 1]) - m0 * S2(w) * D;
    fprintf("u_Y : %.2f %.2f %.2f\n", u_Y(1), u_Y(2), u_Y(3));

    w3d_desired = w3d_des(:, i);
    w3d_desired = - k_w * (w(3) - w3_des(:,i));
    wd_des = [(- nu_(2) / m0 / D(3) + D(1) / D(3) * w3d_desired) ; 
              nu_(1) / m0 / D(3); 
              w3d_desired];
    lambda = nu_(3) + m0 *D(1) * wd_des(2);

    w_int_lim = 10;
    int_wd_des = min(max(int_wd_des + wd_des*dt_sim, -w_int_lim), w_int_lim);
    tau_ = S(w) * (I0 * w) + I0 * (wd_des - alpha *(w - int_wd_des)) - beta * w - tau_hat;

    fprintf("Computed tau_; lambda : %.2f %.2f %.2f; %.2f\n", tau_(1), tau_(2), tau_(3), lambda )
    thrusts = A_eta\[tau_; 0; 0; lambda];
    %thrusts = inv(A)*[tau_; lambda];
    thrusts_hist = [thrusts_hist, thrusts];

    [Ydd, wd_fd, Rd] = forwardDynamics(Y, Yd, w, wd, R, D, f_e, tau_e, thrusts, A_eta, params);
    Y = Y + Yd * dt_sim + 0.5 * Ydd * dt_sim^2 ;
    Yd = Yd + Ydd * dt_sim;
    w = w + wd_fd * dt_sim;
    wd = wd_fd;
    R = R + Rd *dt_sim;
    [U,~,V] = svd(R);
    R = U*V';   
    
    Y_hist = [Y_hist, Y];
    Yd_hist = [Yd_hist, Yd];
    w_hist = [w_hist, wd];
    wd_hist = [wd_hist, wd];
    wd_des_hist = [wd_des_hist, wd_des]; 
    R_hist{i} = R;
    times = [times; i * dt_sim];
end

%% State plot
figure('Position',[100 50 800 900]);
colors = lines(4);

subplot(4,2,1)
hold on
for j = 1:3
    plot(times, Y_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    plot(times, Y_des(j, 1:i), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
end
legend({'$Y_x$', '$Y_x^{\mathrm{des}}$', ...
        '$Y_y$', '$Y_y^{\mathrm{des}}$', ...
        '$Y_z$', '$Y_z^{\mathrm{des}}$'}, ...
        'Interpreter','latex', 'Location', 'best','FontSize', 12);
title('$\mathbf{Y}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
grid on
% 2. Yd plot
subplot(4,2,2)
hold on
for j = 1:3
    plot(times, Yd_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    plot(times, Yd_des(j, 1:i), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
end
legend({'$\dot{Y}_x$', '$\dot{Y}_x^{\mathrm{des}}$', ...
        '$\dot{Y}_y$', '$\dot{Y}_y^{\mathrm{des}}$', ...
        '$\dot{Y}_z$', '$\dot{Y}_z^{\mathrm{des}}$'}, ...
        'Interpreter','latex', 'Location', 'best','FontSize', 12);
title('$\dot{\mathbf{Y}}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
grid on
% 3. w plot
subplot(4,2,3)
hold on
for j = 1:3
    plot(times, w_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
plot(times, w3_des(1, 1:i), '--', 'Color', colors(3,:), 'LineWidth', 2.5);  % w_z des
legend({'$\omega_x$', '$\omega_y$', '$\omega_z$', '$\omega_z^{\mathrm{des}}$'}, ...
       'Interpreter', 'latex', 'Location', 'best','FontSize', 12);
title('${\omega}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
grid on
% 4. wd plot
subplot(4,2,4)
hold on
for j = 1:3
    plot(times, wd_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    plot(times, wd_des_hist(j, 1:i), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
end
legend({'$\dot{\omega}_x$', '$\dot{\omega}_x^{\mathrm{des}}$', ...
        '$\dot{\omega}_y$', '$\dot{\omega}_y^{\mathrm{des}}$', ...
        '$\dot{\omega}_z$', '$\dot{\omega}_z^{\mathrm{des}}$'}, ...
        'Interpreter','latex', 'Location', 'best','FontSize', 12);
title('$\dot{{\omega}}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
grid on

% 5. f plot
subplot(4,2,5)
hold on
for j = 1:3
    plot(times, f_e_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    plot(times, f_hat_hist(j, 1:i), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
end
legend({'$f_x$', '$\hat{f}_x$', ...
        '$f_y$', '$\hat{f}_y$',...
        '$f_z$', '$\hat{f}_z$',}, ...
        'Interpreter','latex', 'Location', 'best','FontSize', 12);
title('$\dot{{\omega}}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
grid on

% 6. f plot
subplot(4,2,6)
hold on
for j = 1:3
    plot(times, tau_e_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    plot(times, tau_hat_hist(j, 1:i), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
end
legend({'$\tau_x$', '$\hat{\tau}_x$', ...
        '$\tau_y$', '$\hat{\tau}_y$',...
        '$\tau_z$', '$\hat{\tau}_z$',}, ...
        'Interpreter','latex', 'Location', 'best','FontSize', 12);
title('$\dot{{\omega}}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
grid on
% 7. thrusts
subplot(4,2,7)
hold on
for j = 1:4
    plot(times, thrusts_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$f_1$', '$f_2$', '$f_3$', '$f_4$'}, ...
        'Interpreter','latex', 'Location', 'best','FontSize', 12);
title('Thrusts', 'Interpreter', 'latex','FontSize', 14)
grid on
%% Error plot
N_plot  = size(Y_hist, 2); 
Y_err   = Y_hist  - Y_des(:, 1:N_plot);
Yd_err  = Yd_hist - Yd_des(:, 1:N_plot);

w3_err   = w_hist(3, :) - w3_des;
wd_err  = wd_hist - wd_des_hist;

f_err = f_e_hist - f_hat_hist;
tau_err = tau_e_hist - tau_hat_hist;

figure('Position',[900 100 800 800]);
colors = lines(3);

subplot(3,2,1);  hold on;
for j = 1:3
    plot(times, Y_err(j,:), 'Color', colors(j,:), 'LineWidth', 1.3);
end
legend({'$e_{Y_x}$','$e_{Y_y}$','$e_{Y_z}$'}, ...
       'Interpreter','latex','Location','best','FontSize', 12);
ylabel('[m]','Interpreter','latex');
title('Position Error $\mathbf e_Y$', 'Interpreter','latex','FontSize', 14);
grid on

subplot(3,2,2);  hold on;
for j = 1:3
    plot(times, Yd_err(j,:), 'Color', colors(j,:), 'LineWidth', 1.3);
end
legend({'$e_{\dot Y_x}$','$e_{\dot Y_y}$','$e_{\dot Y_z}$'}, ...
       'Interpreter','latex','Location','best','FontSize', 12);
ylabel('[$\mathrm{m/s}$]','Interpreter','latex');
title('Velocity Error $\mathbf e_{\dot Y}$', 'Interpreter','latex','FontSize', 14);
grid on

subplot(3,2,3);  hold on;
plot(times, w3_err, 'Color', colors(j,:), 'LineWidth', 1.3);
legend({'$e_{\omega_3}$'}, ...
       'Interpreter','latex','Location','best','FontSize', 12);
ylabel('[$\mathrm{rad/s}$]','Interpreter','latex');
title('Angular-Rate Error $\mathbf e_{\omega}$', 'Interpreter','latex','FontSize', 14);
grid on
% 4) 각가속도 오차
subplot(3,2,4);  hold on;
for j = 1:3
    plot(times, wd_err(j,:), 'Color', colors(j,:), 'LineWidth', 1.3);
end
legend({'$e_{\dot\omega_x}$','$e_{\dot\omega_y}$','$e_{\dot\omega_z}$'}, ...
       'Interpreter','latex','Location','best','FontSize', 12);
ylabel('[$\mathrm{rad/s^2}$]','Interpreter','latex');
title('Angular-Acceleration Error $\mathbf e_{\dot\omega}$', ...
      'Interpreter','latex','FontSize', 14);
grid on
sgtitle('Tracking Errors','Interpreter','latex','FontSize', 16);

% force estimation error
subplot(3,2,5);  hold on;
for j = 1:3
    plot(times, f_err(j,:), 'Color', colors(j,:), 'LineWidth', 1.3);
end
legend({'$e_{f_x}$','$e_{f_y}$','$e_{f_z}$'}, ...
       'Interpreter','latex','Location','best','FontSize', 12);
ylabel('[$\mathrm{N}$]','Interpreter','latex');
title('Force estimation Error $\mathbf e_{f}$', ...
      'Interpreter','latex','FontSize', 14);
grid on
% moment estimation error
subplot(3,2,6);  hold on;
for j = 1:3
    plot(times, tau_err(j,:), 'Color', colors(j,:), 'LineWidth', 1.3);
end
legend({'$e_{\tau_x}$','$e_{\tau_y}$','$e_{\tau_z}$'}, ...
       'Interpreter','latex','Location','best','FontSize', 14);
ylabel('[$\mathrm{N}$]','Interpreter','latex');
title('Moment estimation Error $\mathbf e_{f}$', ...
      'Interpreter','latex','FontSize', 14);
grid on

sgtitle('Tracking Errors','Interpreter','latex','FontSize', 16);
%%
function [Ydd, wd_fd, Rd] = forwardDynamics(Y, Yd, w, wd, R, D, f_e, tau_e,thrusts, A_eta, params)
    m0 = params{1}; I0 = params{2}; mu = params{3}; r= params{4}; gravity = params{16};

    wrench = A_eta * thrusts;
    %A = [r r -r -r;  -r r r -r;mu -mu mu -mu; 0 0 0 0; 0 0 0 0; 1 1 1 1];
    %wrench = A * thrusts;
    tau = wrench(1:3);
    force = wrench(4:6);

    fprintf("Real tau; force : %.2f %.2f %.2f; %.2f %.2f %.2f\n\n", tau(1), tau(2), tau(3), force(1),  force(2),  force(3));

    Ydd = (m0 * R * (S(wd) + S2(w)) * D + R * force - m0 * norm(gravity) * [0; 0; 1] + f_e) / m0;
    wd_fd = inv(I0) * (tau + tau_e - S(w) * I0 * w);
    Rd = R * S(w);
end

function [Y_des, Yd_des, Ydd_des, w3_des, w3d_des] = get_traj(Y_hover, w3_hover, N_sim, dt_sim)
    Y_des = []; Yd_des = []; Ydd_des = []; 
    w3_des = []; w3d_des = [];
    Y_prev = Y_hover; Yd_prev = zeros(3,1);
    w3_prev = 0 / 180 * pi;

    for i = 1:N_sim
        Y = Y_hover;
        Yd = (Y - Y_prev) / dt_sim; 
        Ydd = (Yd - Yd_prev) / dt_sim; 
        Y_des = [Y_des, Y];
        Yd_des = [Yd_des, Yd];
        Ydd_des = [Ydd_des, Ydd];
       
        w3 = w3_hover;
        w3d = (w3 - w3_prev) / dt_sim;
        w3_des = [w3_des, w3];
        w3d_des = [w3d_des, w3d];

        Y_prev = Y; Yd_prev = Yd;
        w3_prev = w3;
    end
end