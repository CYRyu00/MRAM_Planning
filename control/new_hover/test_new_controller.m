addpath("../dynamics", "../functions", "../../params" )
clear; close all
params = define_params_ver2();
%mq = params{1}; Iq = params{2}; 
mu = params{3}; r = params{4}; d = params{5};
thrust_limit= params{6}; gravity = params{16};
l1 = 0.35; l2 = 0.30;
m0 = 2.0;

a = 0.2; b = 0.2;c = 0.1;
J_quad = diag([b^2+c^2, a^2+c^2, a^2+b^2]) /5 * m0 * 0.7;

a = (l1+l2)/2; b = 0.05; c = 0.05;
J_tool = diag([b^2+c^2, a^2+c^2, a^2+b^2]) /5 * m0 * 0.3;

B = [r r -r -r;  -r r r -r; mu -mu mu -mu; 1 1 1 1];

e_1 = [1; 0; 0]; e_2 = [0; 1; 0]; e_3 = [0; 0; 1];
%% inertia
num_AMs = 5;
AM_mass = 0; % mass of shape
AM_inertia = zeros(3, 3); % inertia w.r.t. its com
AM_com = [0; 0; 0];% 0 to com
r_0j = cell(num_AMs, 1); % 0 to j'th module
r_cj = cell(num_AMs, 1); % com to j'th module
I_cj = cell(num_AMs, 1); % inertia of j'th module w.r.t. com of shpe
mass_ams = m0 * ones(num_AMs, 1);
R_shape = cell(1, num_AMs);
shape_pitch =[0 -10 -20 -10 0 0 10];
% shape_pitch =[0 -10 -10 10 -10 -10 10];
% shape_pitch =[0 0 0 0 0 0 0];

for j = 1:length(shape_pitch)
    R_shape{j} = Ry(shape_pitch(j) / 180 *pi);
end


AM_mass = sum(mass_ams);
for j = 1:num_AMs
    if j == 1
        r_0j{j} = [0; 0; 0];
    else
        r_0j{j} = r_0j{j-1} -l2 * R_shape{j-1} * e_1 - l1 * R_shape{j} * e_1;% p_core to j
    end
    AM_com = AM_com + r_0j{j} * mass_ams(j)/AM_mass;
end

% compute AM_inertia
for j = 1:num_AMs
    r_cj{j} = r_0j{j} - AM_com;% p_com to j
    I_cj{j} = J_tool + J_quad + mass_ams(j) * (r_cj{j}' * r_cj{j} * eye(3,3) - r_cj{j} * r_cj{j}'); % TODO: compute Ib It for 3D 
    AM_inertia = AM_inertia + I_cj{j};
end

M = zeros(num_AMs * 6, num_AMs * 6);
for j =1:num_AMs
    M(3*j-2: 3*j, 3*j-2: 3*j) = mass_ams(j) * eye(3);
    M(3*(j+num_AMs)-2: 3*(j+num_AMs), 3*(j+num_AMs)-2: 3*(j+num_AMs)) = J_tool + J_quad;
end
%%
wn = 2; damp = 1.1; % 2, 1.1
kp_M = wn^2; 
kv_M = 2 * damp *sqrt(kp_M);

wn = 2; damp = 1.1; % 1, 1.1
kp_z = wn^2; 
kv_z = 2 * damp *sqrt(kp_z);

kp_M = diag([kp_M, kp_M, kp_z]);
kv_M = diag([kv_M, kv_M, kv_z]);

kp = AM_mass * kp_M;
kv = AM_mass * kv_M;
        

kw_I = 10; % 10 or 20 or 100
kR_I = 10;
kw = kw_I *AM_inertia;
kR = kR_I *AM_inertia;
c2 = kw / 10;
beta2 = 5e1;

% Backstepping gain
epsilon = kv_M(1, 1) * 0.7; % kv_M(1, 1) * 0.7
alpha = 10; % 4 or 10 * num_AMs
gamma = alpha * 1; % alpha * 1.0 
beta = 10.0; % 1~50

% servo moter
ki_servo = 1e0;
kp_servo = 1e2; % 1e1
kd_servo = 1e2; % 1e1
damp_servo = 1e-1; % 1e-1 ~ 1e0

% Simulation Parmeters
dt_sim = 0.0001;
N_sim = 50000;
show_video = true;
save_video = false;
video_speed = 2.0;

% Thrust limit and w_des limit
thrust_limit = thrust_limit * 1.0; % 1 ~ 3
w_des_limit = 2.0; % 2 ~ 1

% disturbance, payload at EE
payload = 3;
rising_time = 3;
mean = 0.0; max_val = 500.0;
sigma = [1.0; 0; 1.0; 0; 0.5; 0] * payload * 9.81 * 0.0;
noise = zeros(6, 1);
disturb_sim = zeros(N_sim, 6);

% Modeling error
mass_uncertainty = 1.00; 
inertia_uncertainty = 1.00;
kinematic_error = 1.1; %TODO

% X, w_estimation error
X_error = zeros(3, num_AMs);
w_error = zeros(3, num_AMs);
sigma_X = 0 / 100; max_X = 0.05;
sigma_w = 0 / 100; max_w = 0.1; 

% Delay
delay_bs = 0.004 / dt_sim;
delay_quad = 0.004 / dt_sim;
delay_cen = 0.01 / dt_sim; % 0.01

rng('shuffle')

%% trajectory
X_hover = [1; 0; 3] * 1e-1;
rpy_hover = [0, 10, 0] / 180 * pi; velocity = [-0.0; 0; 0]; maximum_X = [10.0; 0; 10.0];
[X_des, Xd_des, Xdd_des, Xddd_des, R_e_des, w_e_des, wd_e_des] = get_traj_hover_manip(X_hover, rpy_hover, velocity, maximum_X, N_sim, dt_sim);
% helix
radius = 0.3;  v_z = 0.1;
omega = 2 * pi * 0.3; 
rpyd  = [0.00; -0.00; 0.0] * 2 * pi;
X_hover = [0.1; 0; 0.3]; rpy_hover = [0, 0, 0] / 180 * pi; 
% [X_des, Xd_des, Xdd_des, Xddd_des, R_e_des, w_e_des, wd_e_des] = get_traj_helix_manip_2d(radius, omega, v_z, rpyd, X_hover, rpy_hover, N_sim, dt_sim);
%% Simulation
X = [0; 0; 0]; Xd = [0; 0; 0]; Xdd = [0; 0; 0];
w = [0; 0; 0]; wd = [0; 0; 0];

phi = zeros(num_AMs, 1); phid = zeros(num_AMs, 1); % rad
theta = zeros(num_AMs, 1); thetad = zeros(num_AMs, 1); thetadd = zeros(num_AMs, 1);
R = eye(3, 3);
Rd = R * S(w);

theta_ref = zeros(num_AMs, 1); thetad_ref = zeros(num_AMs, 1);
e_thetai = zeros(num_AMs, 1); thetad_ref_filtered = zeros(num_AMs, 1);
for j = 1:num_AMs
    R_e_j = R * R_shape{j};
    theta_ref(j) = wrapToPi(atan2(R_e_j(1,3), R_e_j(1,1)) - phi(j));
end

lambda = mass_ams * norm(gravity) * 1.0; lambdad = zeros(num_AMs, 1);
w_x_q_des = zeros(num_AMs, 1); w_y_q_des = zeros(num_AMs, 1); wd_q_des = zeros(3, num_AMs);
f_ext_hat = zeros(3, 1);
f_tot_cen = AM_mass * norm(gravity) * e_3;
tau_ext_hat = zeros(3, 1);
tau_tot_cen = zeros(3, 1);

X_hist = []; Xd_hist = []; w_hist = []; wd_hist = [];
% w_des_hist = []; wd_des_hist = []; 
thrusts_hist = []; lambda_hist = [];
e_p_hist = []; e_d_hist = []; e_w_hist = []; e_R_hist = [];
pitch_hist = []; pitch_des_hist = []; phi_hist = [];
nu_e_hist = []; delta_hat_hist = []; delta_tilde_hist = [];

e_theta_hist = []; e_thetad_hist = []; e_pitch_hist = []; 
thetad_ref_hist = []; thetad_hist = [];
disturb_hist = []; f_ext_hat_hist = []; tau_ext_hat_hist = [];
times = [];
R_hist = cell(N_sim, 1);

tic
for i = 1:N_sim
    if mod(i, N_sim/100) == 0
        fprintf("\n%d / %d", i, N_sim)
    end
    tau_tot = zeros(3, 1);%real
    force_tot = zeros(3, 1); %real
    thrusts = [];
    % if mod(i, delay_quad) == 1 || delay_quad == 1
    %     tau = zeros(1, num_AMs);
    % end

    % estimation error
    X_error_dot = randn(3, num_AMs) * sigma_X;
    X_error = X_error + X_error_dot * dt_sim;
    X_error = min(max(X_error, - max_X), max_X);
    w_error_dot = randn(3, num_AMs) * sigma_w;
    w_error = w_error + w_error_dot * dt_sim;
    w_error = min(max(w_error, - max_w), max_w);
    
    % Central controller
    if mod(i, delay_cen) == 1 || delay_cen == 1
        % position
        f_tot_cen = zeros(3, 1);
        for j = 1:num_AMs
            f_tot_cen = f_tot_cen + lambda(j) * Ry(phi(j)) * e_3; 
        end
        e_p  = (1 - X_error(:, j)).*X - X_des(:, i); 
        e_pd = (1 - X_error(:, j)).*Xd - Xd_des(:, i);
        e_pdd_hat = -9.81 * e_3 + (f_tot_cen + f_ext_hat) / AM_mass - Xdd_des(:, i);
        % e_pdd_hat = Xdd - Xdd_des(:, i);
            
        f_tot_des = AM_mass * Xdd_des(:, i) - kv * e_pd - kp * e_p + AM_mass * 9.81 *e_3 - f_ext_hat;
        nu_e = f_tot_cen - f_tot_des;

        f_ext_hatd = beta *(e_pd + epsilon * e_p + kv/AM_mass/gamma* nu_e);
        f_ext_hat = f_ext_hat + f_ext_hatd *dt_sim *delay_cen;

        nud = AM_mass * Xddd_des(:, i) - kv * e_pdd_hat - kp * e_pd - f_ext_hatd; % todo f_ext_hat_d
        fd_tot_des = - alpha * nu_e - gamma * (e_pd + epsilon * e_p) + nud;
        
        % rotaion
        e_R = 2 \ vee(R_e_des{i}' * R - R' * R_e_des{i});
        e_w = w - R' * R_e_des{i} * w_e_des(:, i);
        tau_ext_hatd = beta2 *(e_w + c2 * (AM_inertia \ e_R));
        tau_ext_hat = tau_ext_hat + tau_ext_hatd *dt_sim *delay_cen;
        tau_tot_des = S(w) * AM_inertia * w - kR *e_R -kw *e_w - tau_ext_hat ...
                      - AM_inertia *(S(w) * R' * R_e_des{i} * w_e_des(:, i) - R' * R_e_des{i} * wd_e_des(:, i));

        % optimization : TODO
        % x = [ lambda wx wy 1 ...n , tau_bar 1 ...n] 6n x 1\
        k_lambda = 1;
        k_w_optim = 10;
        k_tau_optim = 1;
        
        options = optimoptions('quadprog', ...
            'Display', 'off');
        H = diag([repmat([k_lambda, k_w_optim, k_w_optim], 1, num_AMs), k_tau_optim * ones(1, 3*num_AMs)]);
        f = zeros(6*num_AMs, 1);
        A_ineq = []; b_ineq = [];
        A_eq = []; b_eq = [];
        lb = -inf; ub = inf;

        % force
        A_tmp = zeros(3, 3*num_AMs);
        for j = 1:num_AMs
            A_tmp(:, 3*j-2:3*j) = Ry(phi(j)) * [0 0 -lambda(j); 0 -lambda(j) 0; 1 0 0];
        end
        A_eq = [A_eq; A_tmp, zeros(3, 3*num_AMs)];
        b_eq = [b_eq; fd_tot_des];
        
        % torque
        A_tmp = repmat(eye(3), 1, num_AMs);
        tau_by_force = zeros(3, 1);
        for j = 1:num_AMs
            tau_by_force = tau_by_force + S(r_cj{j}) * R' * lambda(j) * Ry(phi(j)) * e_3;
        end
        A_eq = [A_eq; zeros(3, 3*num_AMs), A_tmp];
        b_eq = [b_eq; (tau_tot_des - tau_by_force)];
        
        % solution
        sol = quadprog(H, f, A_ineq, b_ineq, A_eq, b_eq, [], [], [], options);

        for j = 1:num_AMs
            w_q_des_prev = [w_x_q_des(j), w_y_q_des(j), 0]; % todo
            lambdad(j) = sol(3*j-2);
            w_x_q_des(j) = min(max(sol(3*j-1), -w_des_limit), w_des_limit);
            w_y_q_des(j) = min(max(sol(3*j), -w_des_limit), w_des_limit);
            tau_bar(j) = sol(3*(j+num_AMs)-1); % for 2D

            lambda(j) = lambda(j) + lambdad(j) *dt_sim *delay_cen;
            wd_q_des(:, j) = ([w_x_q_des(j), w_y_q_des(j), 0] - w_q_des_prev)/dt_sim /delay_cen; 
        end
        %debug
        % nu_e = f_tot_cen - f_tot_des
        % e_w
        % e_R
    end
    
    % Input of Modules
    for j = 1:num_AMs
        
        % rotaion control
        if mod(i, delay_quad) == 1 || delay_quad == 1
            wd_hat = wd; % todo, implement noise
            wd_ref = wd_q_des(:, j) - wd_hat;
            tau(j) = tau_bar(j) - J_quad(2,2) * wd_ref(2);
            
            w_hat = w; % todo, implement noise
            w_ref = [w_x_q_des(j); w_y_q_des(j); 0] - w_hat;
            thetad_ref(j) = w_ref(2);
        end
        % saturated thrusts
        thrust_j = inv(B) * [0; tau(j); 0; lambda(j)];
        thrust_j = min(thrust_limit, max(-thrust_limit, thrust_j));
        
        tau(j) = B(2,:) * thrust_j;
        lambda(j) = B(4, :) * thrust_j;
        thrusts = [thrusts; thrust_j];

        % servo motor
        dt_lpf = 0.001;
        thetad_ref_filtred(j) = (dt_lpf * thetad_ref_filtered(j) + dt_sim * thetad_ref(j)) / (dt_lpf + dt_sim);
        R_e_j = R * R_shape{j};
        w_e_j = R_shape{j}' * w;

        theta(j) = wrapToPi(atan2(R_e_j(1,3), R_e_j(1,1)) - phi(j));
        thetad(j) = w_e_j(2) - phid(j);

        theta_ref(j) = theta_ref(j) + thetad_ref(j) * dt_sim;
        e_theta = theta(j) - theta_ref(j);
        e_thetai(j) = e_thetai(j) + e_theta * dt_sim;
        e_thetad = thetad(j) - thetad_ref(j);
        
        tau_theta(j) = - ki_servo * e_thetai(j) - kp_servo * e_theta - kd_servo * e_thetad - damp_servo * thetad(j); % 1-dim

        % dynamics
        force_tot = force_tot + lambda(j) * Ry(phi(j)) * e_3;
        tau_tot = tau_tot + tau_theta(j) * e_2 +  S(r_cj{j}) * R' * lambda(j) * Ry(phi(j)) * e_3;
    end

    % generate disturbance
    seed_number = i;
    rng(seed_number)
    noise_dot = randn(6, 1) .* sigma;
    noise = noise + noise_dot * dt_sim;
    
    X1 = X + R * r_cj{1};
    X1 = X1 + l1 * R * R_shape{1} * e_1;
    disturb = - payload * min(1, i/rising_time * dt_sim) * 9.81 * [e_3; S(X1 - X) * e_3] + noise;
    disturb = min(max(disturb, - max_val), max_val);
    disturb_sim(i,:) = disturb;

    % Dynamics
    for j = 1:num_AMs
        phidd(j) = J_quad(2,2) \ (tau(j) - tau_theta(j));
        
        phi(j) = phi(j) + phid(j) * dt_sim;
        phid(j) = phid(j) + phidd(j) * dt_sim; 
    end

    % force_tot : sigma lambdai *Ri *e3
    Xdd = (AM_mass * mass_uncertainty) \ (force_tot - (AM_mass * mass_uncertainty) * 9.81 * e_3 + disturb(1:3));
    % tau_tot : sigma tau_i + tau_theta(1) + r x lambda
    wd = inv(AM_inertia * inertia_uncertainty - J_quad * num_AMs) ...
    * (tau_tot - S(w) * (AM_inertia * inertia_uncertainty - J_quad * num_AMs) * w + disturb(4:6));
   
    Xd = Xd + Xdd * dt_sim;
    X = X + Xd * dt_sim + 0.5 * Xdd * dt_sim^2 ;
    w = w + wd * dt_sim;
    Rd = R * S(w);
    R = R + Rd * dt_sim;
    [U_, ~, V_] = svd(R);
    R  = U_ * V_';

    % Log
    times = [times; i * dt_sim];
    X_hist  = [X_hist, X];
    Xd_hist = [Xd_hist, Xd];
    w_hist  = [w_hist, w];
    wd_hist = [wd_hist, wd];
    % wd_des_hist = [wd_des_hist, wd_quad_des]; 
    R_hist{i} = R;
    thrusts_hist = [thrusts_hist, thrusts];
    
    pitch = wrapToPi(atan2(R(1,3), R(1,1)));
    pitch_des = wrapToPi(atan2(R_e_des{i}(1,3), R_e_des{i}(1,1)));
    pitch_hist = [pitch_hist, pitch];
    pitch_des_hist = [pitch_des_hist, pitch_des];
    
    e_p_hist = [e_p_hist, e_p]; 
    e_d_hist = [e_d_hist, e_pd]; 
    e_w_hist = [e_w_hist, e_w];
    e_R_hist = [e_R_hist, e_R];
    nu_e_hist = [nu_e_hist, nu_e];
    f_ext_hat_hist = [f_ext_hat_hist, f_ext_hat];
    tau_ext_hat_hist = [tau_ext_hat_hist, tau_ext_hat];
    disturb_hist = [disturb_hist, disturb];
    
    phi_hist = [phi_hist, phi];

    lambda_hist = [lambda_hist, lambda];
    % 
    % e_pitch_hist = [e_pitch_hist, e_pitch];
    e_theta_hist = [e_theta_hist, e_theta];
    e_thetad_hist = [e_thetad_hist, e_thetad]; 
    thetad_hist = [thetad_hist, thetad];
    thetad_ref_hist = [thetad_ref_hist, thetad_ref];

    % 
    % % Internal force log
    % internal_f_opt_hist = [internal_f_opt_hist, internal_f_opt]; 
    % internal_tau_opt_hist = [internal_tau_opt_hist, internal_tau_opt];
    % 
    % F_ext = [f_ext_hat(1:3); zeros(3*num_AMs -3, 1);
    % e_2 * f_ext_hat(4) - S(r_e) * f_ext_hat(1:3); zeros(3*num_AMs -3, 1)];
    % qd = [zeros(3*num_AMs, 1); repmat(w, num_AMs, 1)];
    % Cqd = zeros(6*num_AMs) * qd;
    % 
    % internal_qd = ((A* (M\ A')) \ Ad) * qd;
    % internal_f_real_qd_hist = [internal_f_real_qd_hist, internal_qd(1:3*(num_AMs-1))]; 
    % internal_tau_real_qd_hist = [internal_tau_real_qd_hist, internal_qd(3*(num_AMs-1)+1:end)];
    % 
    % internal_tq = A_dagger * [zeros(3*num_AMs, 1); input_real(3*num_AMs+1:end)];
    % internal_f_real_tq_hist = [internal_f_real_tq_hist, internal_tq(1:3*(num_AMs-1))]; 
    % internal_tau_real_tq_hist = [internal_tau_real_tq_hist,  internal_tq(3*(num_AMs-1)+1:end)];
    % 
    % internal_real = A_dagger * (input_real + F_ext - Cqd) + internal_qd;
    % internal_f_real_hist = [internal_f_real_hist, internal_real(1:3*(num_AMs-1))]; 
    % internal_tau_real_hist = [internal_tau_real_hist, internal_real(3*(num_AMs-1) + 1:end)];
end
toc
%% State plot
figure('Position',[50 350 500 500]);
colors = lines(6);

subplot(2,2,1)
hold on
for j = [1, 3]
    plot(times, X_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    plot(times, X_des(j, 1:i), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
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
    plot(times, Xd_des(j, 1:i), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
end
legend({'$\dot{X}_x$', '$\dot{X}_x^{\mathrm{des}}$', ...
        '$\dot{X}_z$', '$\dot{X}_z^{\mathrm{des}}$'}, ...
        'Interpreter','latex','FontSize', 10);
title('$\dot{\mathbf{X}}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
grid on
% 3. w plot
subplot(2,2,3)
hold on
for j = 2
    plot(times, w_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    plot(times, w_e_des(j, 1:i), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
end
legend({'$\omega_{y,t}$', '$\omega_{y,t}^{\mathrm{des}}$'}, ...
       'Interpreter', 'latex','FontSize', 12);
title('${\omega_{tool}}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
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
    plot(times, e_w_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$e_w$', 'Interpreter', 'latex','FontSize', 14)
%ylim([-1, 1])
grid on

%nu_e
subplot(3,2,4)
hold on
for j = 1:3
    plot(times, nu_e_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$\nu_{e,i}$', 'Interpreter', 'latex','FontSize', 14)
grid on

%delta_hat
subplot(3,2,5)
hold on
for j = 1:3
    plot(times, f_ext_hat_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    plot(times, disturb_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0, 'LineStyle','--');
end    
plot(times, tau_ext_hat_hist(2, :), 'Color', colors(4,:), 'LineWidth', 1.0);
plot(times, disturb_hist(5, :), 'Color', colors(4,:), 'LineWidth', 1.0, 'LineStyle','--');
    
legend({'$\hat{\Delta}_x^p$','$\Delta_x^p$',...
        '$\hat{\Delta}_y^p$','$\Delta_y^p$',...
        '$\hat{\Delta}_z^p$','$\Delta_z^p$',...
        '$\hat{\Delta}_y^\omega$','$\Delta_y^\omega$'},'Interpreter','latex','FontSize', 12);
title('$\hat{\Delta}$', 'Interpreter', 'latex','FontSize', 14)
grid on

%delta_tilde
subplot(3,2,6)
hold on
for j = 1:3
    if j==4 
        plot(times, disturb_hist(5, :) - f_ext_hat_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    else
        plot(times, disturb_hist(j, :) - f_ext_hat_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    end
end
legend({'$\tilde{\Delta}_x$',...
        '$\tilde{\Delta}_y$',...
        '$\tilde{\Delta}_z$',...
        '$\tilde{\Delta}_y$'},'Interpreter','latex','FontSize', 12);
title('$\tilde{\Delta}$', 'Interpreter', 'latex','FontSize', 14)
grid on

%% Thrusts, servo
figure('Position',[1350 300 400 600]);

% thrusts
subplot(4,1,1)
hold on
plot(times, thrusts_hist, 'LineWidth', 1.0);
ylabel("[N]")
ylim([ 0, thrust_limit*1] )
title('Thrusts', 'Interpreter', 'latex','FontSize', 14)
grid on

subplot(4,1,2)
hold on
plot(times, e_theta_hist / pi * 180, 'LineWidth', 1.0);
ylabel("[degree]")
title('$e\_\theta$', 'Interpreter', 'latex','FontSize', 14)
grid on


subplot(4,1,3)
hold on
plot(times, e_thetad_hist / pi * 180, 'LineWidth', 1.0);
ylabel("[degree]")
title('$e\_\dot\theta$', 'Interpreter', 'latex','FontSize', 14)
grid on

subplot(4,1,4)
hold on
plot(times, thetad_hist / pi * 180, 'LineWidth', 1.0);
plot(times, thetad_ref_hist / pi * 180, 'LineWidth', 1.0 ,'LineStyle', '--');
ylabel("[degree]")
title('$\dot\theta\_ref$', 'Interpreter', 'latex','FontSize', 14)
grid on
%% Video
if show_video
figure('Position',[500 100 500 500]);
dN = 0.1 / dt_sim;
framesPerSecond = 1/dt_sim/dN * video_speed;
rate = rateControl(framesPerSecond);
arrow_len = 0.2;

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
    ylim([X_des(3, i) - (l1+l2)/2*num_AMs*blank, X_des(3, i) + (l1+l2)/2*num_AMs*blank ])
    
    hold on;
    R = R_hist{i};
    
    for j = 1:num_AMs
        R_e_j = R * R_shape{j};
        X_quad = X_hist(:, i) + R * r_cj{j};
        R_quad = Ry(phi_hist(j, i));
        dx = lambda_hist(j, i) / 9.81 * arrow_len * R_quad(1,3);
        dz = (lambda_hist(j, i) - 0 )/ 9.81 * arrow_len * R_quad(3,3);       
        plot(X_quad(1), X_quad(3), 'o', 'Color', 'b');
        quiver(X_quad(1), X_quad(3), dx, dz, 0, 'r', 'LineWidth', 1.5);
        
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
    dx = disturb_hist(1, i) / 9.81 * arrow_len;
    dz = disturb_hist(3, i) / 9.81 * arrow_len;   
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