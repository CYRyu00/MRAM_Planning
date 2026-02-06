addpath("../dynamics", "../functions", "../../params" )
clear; close all
params = define_params_ver2();
mu = params{3}; r = params{4}; d = params{5};
thrust_limit = params{6}; gravity = params{16};
mb = params{17}; cb = params{18}; Ib = params{19};
mt = params{20}; ct = params{21}; It = params{22};
ma = params{23}; ca = params{24}; Ia = params{25};
m0 = mt + ma + mb;

B = [1 1 1 1; r r -r -r;  -r r r -r; mu -mu mu -mu];
l1 = 0.35; l2 = 0.35; % ca = 0.24

e_1 = [1; 0; 0]; e_2 = [0; 1; 0]; e_3 = [0; 0; 1];

a = 0.2; b_ineq = 0.2;c = 0.1;
J_quad = diag([b_ineq^2+c^2, a^2+c^2, a^2+b_ineq^2]) /5 * m0 * 0.7;

a = (l1+l2)/2; b_ineq = 0.05; c = 0.05;
J_tool = diag([b_ineq^2+c^2, a^2+c^2, a^2+b_ineq^2]) /5 * m0 * 0.3;
%% Compute Inertia
num_AMs = 5;
J_c = zeros(3, 3); % inertia w.r.t. its com
AM_com = [0; 0; 0];% 0 to com
r_0j = cell(num_AMs, 1); % 0 to j'th module
r_cj = cell(num_AMs, 1); % com to j'th module
I_cj = cell(num_AMs, 1); % inertia of j'th module w.r.t. com of shpe
mass_ams = m0 * ones(num_AMs, 1);
R_shape = cell(1, num_AMs);
shape_pitch = [0 -10 -10 10 -10 -10 10]; %[0 -10 -10 10 -10 -10 10]

for j = 1:length(shape_pitch)
    R_shape{j} = Ry(shape_pitch(j) / 180 *pi);
end

m_c = sum(mass_ams);
for j = 1:num_AMs
    if j == 1
        r_0j{j} = [0; 0; 0];
    else
        r_0j{j} = r_0j{j-1} -l2 * R_shape{j-1} * e_1 - l1 * R_shape{j} * e_1;% p_core to j
    end
    AM_com = AM_com + r_0j{j} * mass_ams(j)/m_c;
end

% compute AM_inertia
for j = 1:num_AMs
    r_cj{j} = r_0j{j} - AM_com;% p_com to j
    I_cj{j} = J_tool + J_quad + mass_ams(j) * (r_cj{j}' * r_cj{j} * eye(3,3) - r_cj{j} * r_cj{j}'); % TODO: compute J_tool, J_quad for 3D 
    J_c = J_c + I_cj{j};
end
J_tool_tot = J_c - num_AMs * J_quad;

%% Set Parameters
% control gains
wn = 1; damp = 1.1; % 2, 1.1
kp_M = wn^2; 
kv_M = 2 * damp *sqrt(kp_M);
wn = 1; damp = 1.1; % 2, 1.1
kp_z = wn^2; 
kv_z = 2 * damp *sqrt(kp_z);
k_p = diag([kp_M, kp_M, kp_z]) * m_c;
k_v = diag([kv_M, kv_M, kv_z]) * m_c;

wn = 1; damp = 1.5; % 2, 1.1
kR_M = wn^2; 
kw_M = 2 * damp *sqrt(kR_M);
k_R = diag([kR_M, kR_M, kR_M]) * J_c;
k_w = diag([kw_M, kw_M, kw_M]) * J_c;

beta_t = 10.0; % 10 ~ 100
epsilon_t = 0.2 * k_v / m_c; % 0.2
beta_r = 10.0; % 10 ~ 100
epsilon_r = 0.2 * k_w / J_c;

% servo moter
kp_servo = 0.4 * 1 * 180 / pi; %  0.4 * 5 * 180/pi : 1.0Nm per 2.5 degree  [Nm/rad] 
kd_servo = kp_servo^0.5 * 1.2; % kp_servo^0.5 * 1.2 
damp_servo = kd_servo * 0.2; % kd_servo * 0.2 

% saturation parmeters
thrust_limit = thrust_limit * 3.0; % 1 ~ 3
w_des_limit = 2.0; % 2 ~ 1
thetad_limit = 1.0;

% central controller
options = optimoptions('quadprog', ...
    'Algorithm', 'interior-point-convex', ...
    'Display', 'off', ...
    'OptimalityTolerance', 1e-5, ...
    'MaxIterations', 100, ...
    'ConstraintTolerance', 1e-5);
lb = -thrust_limit * ones(4*num_AMs, 1);
ub = thrust_limit * ones(4*num_AMs, 1);

% Simulation Parmeters
dt_sim = 0.001;
N_sim = 10000;
show_video = true;
save_video = false;
video_speed = 1.0;

% Control freq.
step_central_control = 1 /dt_sim /25; % 1/ dt_sim/ Hz
step_config_opt = 1 /dt_sim /10;
delay_communication = 0.02 / dt_sim; % 20ms

% Disturbance, payload at EE
payload = 0.0;
rising_time = 3;
mean = 0.0; max_val = 500.0;
sigma = [1.0; 0; 1.0; 0; 0.5; 0] * payload * 9.81 * 0.0;
noise = zeros(6, 1);
disturb_hist = zeros(N_sim, 6);

% Modeling error
mass_uncertainty = 1.00; 
inertia_uncertainty = 1.00;

% X, w_estimation error
X_error = zeros(3, num_AMs);
w_error = zeros(3, num_AMs);
sigma_X = 0 / 100; max_X = 0.05;
sigma_w = 0 / 100; max_w = 0.1; 
rng('shuffle')

%% Trajectory
% linear
X_hover = [1; 0; 3] * 1e-1;
rpy_hover = [0, 10, 0] / 180 * pi; velocity = [-0.1; 0; 0]; maximum_X = [10.0; 0; 10.0];
[X_des, Xd_des, Xdd_des, Xddd_des, R_e_des, w_e_des, wd_e_des] = get_traj_hover_manip(X_hover, rpy_hover, velocity, maximum_X, N_sim, dt_sim);
% helix
radius = 0.3;  v_z = 0.0;
omega = 2 * pi * 0.2; 
rpyd  = [0.00; -0.00; 0.0] * 2 * pi;
X_hover = [0.1; 0; 0.3]; rpy_hover = [0, 10, 0] / 180 * pi; 
[X_des, Xd_des, Xdd_des, Xddd_des, R_e_des, w_e_des, wd_e_des] = get_traj_helix_manip_2d(radius, omega, v_z, rpyd, X_hover, rpy_hover, N_sim, dt_sim);
%% Simulation
X = [0; 0; 0]; Xd = [0; 0; 0]; Xdd = [0; 0; 0];
w = [0; 0; 0]; wd = [0; 0; 0];
R = eye(3, 3); Rd = R * S(w);

phi_q = (2*rand(num_AMs, 1) - 1) * 10 / 180 *pi; phid_q = zeros(num_AMs, 1);
R_q = cell(1, num_AMs);
for j = 1:num_AMs
    R_q{j} = Ry(phi_q(j));
end
theta = zeros(num_AMs, 1); thetad = zeros(num_AMs, 1); thetadd = zeros(num_AMs, 1);
theta_ref = zeros(N_sim, num_AMs);
X_hist = []; Xd_hist = []; R_hist = cell(N_sim, 1); w_hist = []; wd_hist = [];
pitch_hist = []; pitch_des_hist = []; phi_q_hist = [];
e_p_hist = []; e_d_hist = []; e_R_hist = []; e_w_hist = [];
theta_hist = []; theta_ref_hist = []; theta_opt_hist = [];
e_thetaj_hist = []; e_thetadj_hist = []; 
thrusts_hist = [];
wrench_ext_tilde_hist = [];
times = [];

tic
for i = 1:N_sim
    if mod(i, N_sim/100) == 0
        fprintf("\n%d / %d\n", i, N_sim)
    end
    for j = 1:num_AMs
        R_e_j = R * R_shape{j}; % j'th module's tool rotation
        w_e_j = R_shape{j}' * w;
        thetad(j) = wrapToPi(w_e_j(2) - phid_q(j));
        theta(j) = wrapToPi(atan2(R_e_j(1,3), R_e_j(1,1)) - phi_q(j));
    end
    if i == 1
        theta_ref(1, :) = theta;
        wrench_ext_hat = zeros(6, 1); % [force; moment]
    end
    % estimation error
    X_error_dot = randn(3, num_AMs) * sigma_X;
    X_error = X_error + X_error_dot * dt_sim;
    X_error = min(max(X_error, - max_X), max_X);
    w_error_dot = randn(3, num_AMs) * sigma_w;
    w_error = w_error + w_error_dot * dt_sim;
    w_error = min(max(w_error, - max_w), max_w);
    
    % Configuration optimization : TODO
    if mod(i, step_config_opt) == 1 || step_config_opt == 1
        theta_opt = ([-2 -1 0 1 2]' * -30 + shape_pitch(1:num_AMs)' + 0 * sin(2 * pi * 0.2 * i * dt_sim)* [1 1 1 1 1]') / 180 *pi;
    end

    % Central controller
    if mod(i, step_central_control) == 1 || step_central_control == 1 
        e_p  = (1 - X_error(:, 1)).*X - X_des(:, i); 
        e_pd = (1 - X_error(:, 1)).*Xd - Xd_des(:, i);
        force_ext_hatd = beta_t * (e_pd + epsilon_t * e_p);
        wrench_ext_hat(1:3) = wrench_ext_hat(1:3) + force_ext_hatd * dt_sim * step_central_control; 
        force_des = m_c * Xdd_des(:, i) - k_p * e_p - k_v * e_pd + m_c * 9.81 * e_3 - wrench_ext_hat(1:3);
        
        R_hat = R;
        e_R = 0.5 * vee(R_e_des{i}' * R_hat - R_hat' * R_e_des{i});
        e_w = (1 - w_error(:, 1)).* w - R_hat' * R_e_des{i} * w_e_des(:, i);
        torque_ext_hatd = beta_r * (e_w + epsilon_r * e_R); %TODO
        wrench_ext_hat(4:6) = wrench_ext_hat(4:6) + torque_ext_hatd * dt_sim * step_central_control; 
        tau_des = J_c * wd_e_des(:, i) - k_R * e_R - k_w * e_w - wrench_ext_hat(4:6);
        
        % Thrusts Optimization - QP
        H = diag(ones(4*num_AMs, 1)); % 0.5 x' * H * x + f' * x
        f = zeros(4*num_AMs, 1); 
        
        A_ineq = [];
        b_ineq = [];
        
        A_force = [];
        A_tau = [];
        for j = 1:num_AMs
            R_tq_j = R' * R_q{j};
            A_force = [A_force, R_tq_j * e_3, zeros(3, 3)];
            A_tau = [A_tau, S(r_cj{j}) * R_tq_j * e_3, eye(3, 3)];
        end
        A_eq = [A_force; A_tau] * kron(eye(num_AMs), B);
        b_eq = [R' * force_des; tau_des];
        
        if i == 1
            sol_init = zeros(4*num_AMs, 1);
        end

        
        [sol, fval, exitflag, output] = quadprog(H, f, A_ineq, b_ineq, A_eq, b_eq, lb, ub, sol_init, options);
        
        sol_init = sol;
        if exitflag ~= 1 && exitflag ~= 2 
            fprintf("step : %d, exitflag %d\n", i, exitflag)
            disp(output)
        end

        % Servo command
        for timestep = 1:step_central_control
            thetad_ref(i + timestep - 1, :) = min(thetad_limit, max(- thetad_limit, ... 
            (theta_opt - theta) / (dt_sim * step_config_opt) )); %TODO :뭔가 array 형식으로 넘겨줘야함
            theta_ref(i + timestep, :) = theta_ref(i + timestep - 1, :) + thetad_ref(i + timestep - 1, :) * dt_sim;    
        end
    end
    if i == 1
        sols = [];
    end

    sols = [sols, sol];
    
    % generate disturbance
    seed_number = i;
    rng(seed_number)
    noise_dot = randn(6, 1) .* sigma;
    noise = noise + noise_dot * dt_sim;
    
    X1 = X + R * r_cj{1} + l1 * R * R_shape{1} * e_1;
    disturb = - payload * min(1, i/rising_time * dt_sim) * 9.81 * [e_3; S(X1 - X) * e_3] + noise;
    disturb = min(max(disturb, - max_val), max_val);
    disturb_hist(i, :) = disturb;

    tau_tot = zeros(3, 1);
    force_tot = zeros(3, 1);
    thrusts = [];

    for j = 1:num_AMs
        if i - delay_communication > 0 
            thrust_j = sols(4*j-3: 4*j, i - delay_communication); %TODO
        else
            thrust_j = zeros(4, 1); %sols(4*j-3: 4*j, 1);
        end
        thrust_j = min(thrust_limit, max(-thrust_limit, thrust_j));
        tau_q(j) = B(3,:) * thrust_j;
        lambda(j) = B(1, :) * thrust_j;
        thrusts = [thrusts; thrust_j];

        % servo motor
        if i - delay_communication > 0 
            theta_ref_servo = theta_ref(i - delay_communication, j);
            thetad_ref_servo = thetad_ref(i - delay_communication, j);
        else
            theta_ref_servo = theta_ref(i, j);
            thetad_ref_servo = thetad_ref(i, j);
        end
        e_theta = theta(j) - theta_ref_servo;
        e_thetad = thetad(j) - thetad_ref_servo;
        tau_theta(j) = - kp_servo * e_theta - kd_servo * e_thetad - damp_servo * thetad(j); % 1-dim
        
        term = S(r_cj{j}) * R' * R_q{j} * [0; 0; lambda(j)]; % 3-dim
        tau_tot(2) = tau_tot(2) + tau_theta(j) + term(2); % 1-dim
        force_tot = force_tot + lambda(j) * R_q{j} * e_3; % 3-dim
    end

    % Dynamics
    % force_tot : sigma lambda(j) * Rqj *e3
    Xdd = (m_c * mass_uncertainty) \ (force_tot - (m_c * mass_uncertainty) * 9.81 * e_3 + disturb(1:3));
    % tau_tot : sigma tau_theta(j) + r x lambda
    wd = inv(J_tool_tot * inertia_uncertainty) * (tau_tot - S(w) * (J_tool_tot * inertia_uncertainty) * w + disturb(4:6));

    Xd = Xd + Xdd * dt_sim;
    X = X + Xd * dt_sim;
    w = w + wd * dt_sim;
    Rd = R * S(w);
    R = R + Rd * dt_sim; [U_, ~, V_] = svd(R); R = U_ * V_';

    for j = 1:num_AMs
        phidd_q = J_quad(2,2) \ (tau_q(j) - tau_theta(j));
        phid_q(j) = phid_q(j) + phidd_q * dt_sim;
        phi_q(j) = phi_q(j) + phid_q(j) * dt_sim;
        R_q{j} = Ry(phi_q(j));
        
        R_e_j = R * R_shape{j}; % j'th module's tool rotation
        w_e_j = R_shape{j}' * w;
        thetad(j) = wrapToPi(w_e_j(2) - phid_q(j));
        theta(j) = wrapToPi(atan2(R_e_j(1,3), R_e_j(1,1)) - phi_q(j));
    end

    % Log
    X_hist  = [X_hist, X];
    Xd_hist = [Xd_hist, Xd];
    w_hist  = [w_hist, w];
    wd_hist = [wd_hist, wd];
    R_hist{i} = R;

    e_p_hist = [e_p_hist, X - X_des(:, i)]; 
    e_d_hist = [e_d_hist, Xd - Xd_des(:, i)]; 
    e_w_hist = [e_w_hist, w - R' * R_e_des{i} * w_e_des(:, i)];
    e_R_hist = [e_R_hist, 0.5 * vee(R_e_des{i}' * R - R' * R_e_des{i})];
    e_thetaj_hist = [e_thetaj_hist, e_theta]; 
    e_thetadj_hist = [e_thetadj_hist, e_thetad];
    
    thrusts_hist = [thrusts_hist, thrusts];

    times = [times, i * dt_sim];
    pitch = atan2(R(1,3), R(1,1));
    pitch_des = atan2(R_e_des{i}(1,3), R_e_des{i}(1,1));
    pitch_hist = [pitch_hist, pitch];
    pitch_des_hist = [pitch_des_hist, pitch_des]; 
    phi_q_hist = [phi_q_hist, phi_q];
    theta_hist = [theta_hist, theta];
    theta_opt_hist = [theta_opt_hist, theta_opt];
    wrench_ext_tilde_hist = [wrench_ext_tilde_hist, disturb - wrench_ext_hat];
end
toc
%% State plot
figure('Position',[50 350 500 500]);
colors = lines(6);

subplot(2,2,1)
hold on
for j = [1, 3]
    plot(times, X_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    plot(times, X_des(j, :), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
end
legend({'$X_x$', '$X_x^{\mathrm{des}}$', ...
        '$X_z$', '$X_z^{\mathrm{des}}$'}, ...
        'Interpreter','latex','FontSize', 10);
title('$\mathbf{X}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
grid on
% 2. Xd plot
subplot(2,2,2)
hold on
for j = [1, 3]
    plot(times, Xd_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    plot(times, Xd_des(j, :), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
end
legend({'$\dot{X}_x$', '$\dot{X}_x^{\mathrm{des}}$', ...
        '$\dot{X}_z$', '$\dot{X}_z^{\mathrm{des}}$'}, ...
        'Interpreter','latex','FontSize', 10);
title('$\dot{\mathbf{X}}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
grid on
% 3. w plot
subplot(2,2,3)
hold on
plot(times, w_hist(2, :), 'Color', colors(j,:), 'LineWidth', 1.0);
plot(times, w_e_des(2, :), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
legend({'$\omega_{y,i}$', '$\omega_{y,i}^{\mathrm{des}}$'}, ...
       'Interpreter', 'latex','FontSize', 12);
title('${\omega_i}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
%ylim([-1, 1])
grid on

% e_pitch
subplot(2,2,4)
hold on
plot(times, pitch_hist / pi * 180, 'Color', colors(2,:), 'LineWidth', 1.0);
plot(times, pitch_des_hist / pi * 180, '--', 'Color', colors(2,:), 'LineWidth', 2.5);
ylabel("degree")
legend({'$\phi_{e,1}$', '$\phi_{e,1}^{\mathrm{des}}$'}, ...
       'Interpreter', 'latex','FontSize', 12);
title('${\phi_{e}}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
grid on
%% Error 
figure('Position',[750 350 600 600]);
colors = lines(6);

%e_p
subplot(3,2,1)
hold on
for j = 1:3
    plot(times, e_p_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_p$', 'Interpreter', 'latex','FontSize', 14)
grid on

%e_d
subplot(3,2,2)
hold on
for j = 1:3
    plot(times, e_d_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_v$', 'Interpreter', 'latex','FontSize', 14)
grid on

%e_w
subplot(3,2,3)
hold on
for j = 1:3
    plot(times, e_R_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_R$', 'Interpreter', 'latex','FontSize', 14)
%ylim([-1, 1])
grid on

%e_w
subplot(3,2,4)
hold on
for j = 1:3
    plot(times, e_w_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_w$', 'Interpreter', 'latex','FontSize', 14)
%ylim([-1, 1])
grid on

% force_tilde
subplot(3,2,5)
hold on
for j = 1:3
    plot(times, wrench_ext_tilde_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$\tilde{f}^{ext}$', 'Interpreter', 'latex','FontSize', 14)
grid on

% tau_tilde
subplot(3,2,6)
hold on
for j = 4:6
    plot(times, wrench_ext_tilde_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$\tilde\tau^{ext}$', 'Interpreter', 'latex','FontSize', 14)
grid on

%% force/disturbance
figure('Position',[1350 550 500 400]);

% thrusts
subplot(4,1,1)
hold on
plot(times, thrusts_hist, 'LineWidth', 1.0);
%ylim([ -thrust_limit, thrust_limit] )
title('Thrusts', 'Interpreter', 'latex','FontSize', 14)
grid on

%% End effector Error 
figure('Position',[1350 50 500 450]);
colors = lines(6);

subplot(3,1,1)
hold on
plot(times, theta_hist(1, :) / pi * 180, 'Color', colors(1,:), 'LineWidth', 1.0);
plot(times, theta_ref(1:N_sim, 1) / pi * 180, 'Color', colors(2,:), 'LineWidth', 1.0, 'LineStyle','--');
plot(times, theta_opt_hist(1, :) / pi * 180, 'Color', colors(3,:), 'LineWidth', 1.0, 'LineStyle',':');
ylabel("[degree]")
legend({'$\theta$','$\theta^{ref}$','$\theta^{opt}$'},'Interpreter','latex','FontSize', 12);
title('$e_\theta$', 'Interpreter', 'latex','FontSize', 14)
grid on

subplot(3,1,2)
hold on
plot(times, e_thetaj_hist(1, :) / pi * 180, 'Color', colors(1,:), 'LineWidth', 1.0);
ylabel("[degree]")
legend({'$pitch$'},'Interpreter','latex','FontSize', 12);
title('$e_\theta$', 'Interpreter', 'latex','FontSize', 14)
grid on

subplot(3,1,3)
hold on
plot(times, e_thetadj_hist(1, :) / pi * 180, 'Color', colors(1,:), 'LineWidth', 1.0);
ylabel("[degree/sec]")
legend({'$pitch$'},'Interpreter','latex','FontSize', 12);
title('$e_{\dot\theta}$', 'Interpreter', 'latex','FontSize', 14)
grid on
%% Video
if show_video
figure('Position',[500 100 700 400]);
dN = 0.1 / dt_sim;
framesPerSecond = 1/dt_sim/dN * video_speed;
rate = rateControl(framesPerSecond);
arrow_len = 0.4;

if save_video
    video_filename = '../images/test.avi';
    video = VideoWriter(video_filename);
    video.FrameRate = framesPerSecond;
    open(video);
end

for i = 1:dN:N_sim
    clf;
    grid on; axis equal;
    set(gca, 'GridLineWidth', 1.0, 'GridAlpha', 0.3);
    blank = 1.5;
    xlim([X_des(1, i) - (l1+l2)/2*num_AMs*blank, X_des(1, i) + (l1+l2)/2*num_AMs*blank ] )
    ylim([X_des(3, i) - (l1+l2)/2*num_AMs*blank/2, X_des(3, i) + (l1+l2)/2*num_AMs*blank/2])
    
    hold on;
    R = R_hist{i};
    
    for j = 1:num_AMs
        R_e_j = R * R_shape{j};
        X_quad = X_hist(:, i) + R * r_cj{j};
        R_quad = Ry(phi_q_hist(j, i));

        tau_q_j = B(3,:) * thrusts_hist(4*j-3:4*j, i);
        lambda_j = B(1, :) * thrusts_hist(4*j-3:4*j, i);

        plot(X_quad(1), X_quad(3), 'o', 'Color', 'b');
        
        X1 = X_quad + r * R_quad * e_1;
        X2 = X_quad - r * R_quad * e_1;
        plot([X1(1), X2(1)], [X1(3), X2(3)], 'b-', 'LineWidth', 2.0);

        dx = thrusts_hist(4*j-3, i) / 9.81 * arrow_len * R_quad(1,3);
        dz = thrusts_hist(4*j-3, i) / 9.81 * arrow_len * R_quad(3,3);       
        quiver(X1(1), X1(3), dx, dz, 0, 'r', 'LineWidth', 1.5);
        dx = thrusts_hist(4*j-2, i) / 9.81 * arrow_len * R_quad(1,3);
        dz = thrusts_hist(4*j-2, i) / 9.81 * arrow_len * R_quad(3,3);       
        quiver(X2(1), X2(3), dx, dz, 0, 'r', 'LineWidth', 1.5);

        X1 = X_quad + l1 * R_e_j * e_1;
        X2 = X_quad - l2 * R_e_j * e_1;
        plot([X1(1), X2(1)], [X1(3), X2(3)], 'b-', 'LineWidth', 1.0);
   
    end

    R_des = R_e_des{i};
    
    for j = 1:num_AMs
        R_e_des_j = R_des * R_shape{j};
        X_des_quad = X_des(:, i) + R_des * r_cj{j};      
        plot(X_des_quad(1), X_des_quad(3), 'o', 'Color', 'black', 'LineWidth', 2.0);

        X1 = X_des_quad + l1 * R_e_des_j * e_1;
        X2 = X_des_quad - l2 * R_e_des_j * e_1;
        plot([X1(1), X2(1)], [X1(3), X2(3)], 'k--', 'LineWidth', 2.0);
    end

    % com
    plot(X_hist(1, i), X_hist(3, i), 'o', 'Color', 'b', 'MarkerSize', 10);
    plot(X_des(1, i), X_des(3, i), 'o', 'Color', 'k', 'MarkerSize', 8);
    
    % payload
    X1 = X_hist(:, i) + R * r_cj{1};
    X1 = X1 + l1 * R * R_shape{1} * e_1;
    dx = disturb_hist(i, 1) / 9.81 * arrow_len;
    dz = disturb_hist(i, 3) / 9.81 * arrow_len;   
    quiver(X1(1), X1(3), dx, dz, 0, 'g', 'LineWidth', 1.5);

    xlabel('X position'); ylabel('Z position');
    title_string = sprintf("time : %.2f sec", times(i));
    title(title_string);
    drawnow

    if save_video
        frame = getframe(gcf);
        writeVideo(video, frame);
        waitfor(rate);
    end
end
if save_video
    close(video)
    fprintf("\n Video saved at %s\n", video_filename);
end
end
%%
function out = Rx(q)
    out = [1, 0, 0;
          0, cos(q), -sin(q);
          0, sin(q), cos(q)];
end
function out = Ry(q)
    out = [cos(q), 0, sin(q);
           0, 1, 0;
           -sin(q), 0, cos(q)];
end