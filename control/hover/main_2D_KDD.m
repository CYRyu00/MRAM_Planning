addpath("../dynamics", "../functions", "../../params" )
clear; close all
params = define_params_ver2();
%mq = params{1}; Iq = params{2}; 
mu = params{3}; r = params{4}; d = params{5};
thrust_limit= params{6}; gravity = params{16};
mb = params{17}; cb = params{18}; Ib = params{19};
mt = params{20}; ct = params{21}; It = params{22};
ma = params{23}; ca = params{24}; Ia = params{25};

% totoal robotic arm
Ib = Ia + Ib;% Ib = Ib * 0.01;
mb = ma + mb;
m0 = mt + mb;

dt_sim = 0.002;
N_sim = 10000;

B = [r r -r -r;  -r r r -r; mu -mu mu -mu; 1 1 1 1];

e_1 = [1; 0; 0];
e_2 = [0; 1; 0];
e_3 = [0; 0; 1];
%% inertia
num_AMs = 4;
AM_mass = 0; % mass of shape
AM_inertia = zeros(3, 3); % inertia w.r.t. its com
AM_com = [0; 0; 0];% 0 to com
r_0j = cell(num_AMs, 1); % 0 to j'th module
r_cj = cell(num_AMs, 1); % com to j'th module
I_cj = cell(num_AMs, 1); % inertia of j'th module w.r.t. com of shpe
mass_ams = m0 * ones(num_AMs, 1);
R_shape = cell(1, num_AMs);
R_shape{1} = eye(3,3);
R_shape{2} = Ry(10 / 180 *pi);
R_shape{3} = Ry(20 / 180 *pi);
R_shape{4} = Ry(30 / 180 *pi);
R_shape{5} = Ry(40 / 180 *pi);
l1 = 0.3; l2 = 0.3; % ca = 0.24

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
    I_cj{j} = Ib + It + mass_ams(j) * (r_cj{j}' * r_cj{j} * eye(3,3) - r_cj{j} * r_cj{j}'); % TODO: compute Ib It for 3D 
    AM_inertia = AM_inertia + I_cj{j};
end

%%
wn = 0.5; damp = 1.2; % 0.5, 1.2
kp_M = wn^2; 
kv_M = 2 * damp *sqrt(kp_M);

wn = 1.0; damp = 1.2; % 1, 1.2
kp_z = wn^2; 
kv_z = 2 * damp *sqrt(kp_z);

kp_M = diag([kp_M, kp_M, kp_z]);
kv_M = diag([kv_M, kv_M, kv_z]);

kw_I = 50; % 10 or 20 or 100

% servo moter
kp_servo = 0.1; % 0.1 / 0.1
kd_servo = 1.0; % 1.0 / 0.1 / 0.03
damp_servo = 0.05; % 0.05

% damped pseudoinverse
damped = 0.0; % 0.3
k_pitch = 0.5; % 1.0

epsilon = kv_M(1, 1) * 0.3; % kv_M(1, 1) * 0.7
alpha = 10; % 4 or 5
gamma = alpha * 1.0; % alpha * 1.0 
N_sim_tmp = 10000;

mass_uncertainty = 1.00; 
inertia_uncertainty = 1.00;
thrust_limit = thrust_limit * 3;

%disturbance
sigma = 0.0; mean = 1.0; max_val = 30.0;
disturb_sim = zeros(N_sim, 6);
disturb = mean * [0; 0.0; 0; 0.0; 0.0; 0.0];
% X, w_estimation error
X_error = zeros(3, num_AMs);
w_error = zeros(3, num_AMs);
sigma_X = 0 / 100; max_X = 0.05;
sigma_w = 0 / 100; max_w = 0.1; 
delay_bs = 0.004 / dt_sim;
delay_quad = 0.004 / dt_sim;

rng('shuffle')

X_hover = [1; 0; 3] * 1e-1;
rpy_hover = [0, 5, 0] / 180 * pi;
[X_des, Xd_des, Xdd_des, Xddd_des, R_e_des, w_e_des, wd_e_des] = get_traj_hover_manip(X_hover, rpy_hover, N_sim, dt_sim);
% helix
radius = 0.2;  v_z = 0.05;
omega = 2 * pi * 0.1; 
rpyd  = [0.00; 0.01; 0.0] * 2 * pi;
X_hover = [0.1; 0; 0.3]; rpy_hover = [0, 0, 0] / 180 * pi; 
[X_des, Xd_des, Xdd_des, Xddd_des, R_e_des, w_e_des, wd_e_des] = get_traj_helix_manip_2d(radius, omega, v_z, rpyd, X_hover, rpy_hover, N_sim, dt_sim);

%%
X = [0; 0; 0]; Xd = [0; 0; 0]; Xdd = [0; 0; 0];
w = [0; 0; 0]; wd = [0; 0; 0];
Xddd_des_quad = zeros(3, num_AMs);

phi = zeros(num_AMs, 1); phid = zeros(num_AMs, 1);
theta = zeros(num_AMs, 1); thetad = zeros(num_AMs, 1); thetadd = zeros(num_AMs, 1);
theta_ref = zeros(num_AMs, 1); %thetad_ref = zeros(2, num_AMs);
R = eye(3, 3);
Rd = R * S(w);
Rt = cell(1, num_AMs);
for j = 1:num_AMs
    Rt{j} = Ry(phi(j));
end

force_prev = mass_ams * norm(gravity) * 1.0;
delta_hat = zeros(3, num_AMs);
w_quad_des_prev = zeros(3, num_AMs);

X_hist = []; Xd_hist = []; w_hist = []; wd_hist = []; R_hist = cell(N_sim, 1); 
w_des_hist = []; wd_des_hist = []; thrusts_hist = []; F_hat_hist = [];
e_p_hist = []; e_d_hist = []; e_w_hist = []; e_psi_hist = [];
nu_e_hist = []; delta_hat_hist = []; delta_tilde_hist = [];
force_per_M_hist = []; delta_hat_x_per_M_hist = [];
e_theta_hist = []; e_thetad_hist = [];
e_pitch_hist = []; e_R_hist = [];
pitch_hist = [];pitch_des_hist = [];
phi_hist = [];
times = [];

for i = 1:N_sim_tmp
    % Semi-decentralized control
    tau_tot = zeros(3, 1);
    tau = zeros(1, num_AMs);
    tau_theta = zeros(num_AMs, 1);
    force_tot = [0; 0; 0];
    thrusts = [];
    force_per_M = [];
    delta_hat_x_per_M = [];
    %estimation error
    X_error_dot = randn(3, num_AMs) * sigma_X;
    X_error = X_error + X_error_dot * dt_sim;
    X_error = min(max(X_error, - max_X), max_X);
    w_error_dot = randn(3, num_AMs) * sigma_w;
    w_error = w_error + w_error_dot * dt_sim;
    w_error = min(max(w_error, - max_w), max_w);

    for j = num_AMs:-1:1
    %for j = 1:num_AMs
        % position control - backstepping
        mj = mass_ams(j);
        rj = r_cj{j};
        kp_j = mj * kp_M;
        kv_j = mj * kv_M;

        R_quad = Rt{j};
        w_quad_j = (1 - w_error(:, j)).*[0; phid(j); 0];
        
        X_quad = X + R * rj; 
        X_des_quad = X_des(:, i) + R_e_des{i} * rj;
        Xd_quad = Xd + Rd * rj; 
        Xd_des_quad = Xd_des(:, i) + R_e_des{i} * S(w_e_des(:, i)) * rj;
        Xdd_des_quad = Xdd_des(:, i) + R_e_des{i}* S2(w_e_des(:, i)) * rj + R_e_des{i}* S(wd_e_des(:,i)) * rj;

        if i < N_sim
            Xddd_des_quad(:, j) = (Xdd_des(:, i+1) + R_e_des{i+1}* S2(w_e_des(:,i)) * rj + R_e_des{i+1}* S(wd_e_des(:,i)) * rj...
                                   - Xdd_des_quad) / dt_sim;
        end
        
        if mod(i, delay_bs) == 1 || delay_bs == 1
        e_p  = (1 - X_error(:, j)).*X_quad - X_des_quad; 
        e_pd = (1 - X_error(:, j)).*Xd_quad - Xd_des_quad;
        
        e_pdd_hat = gravity + force_prev(j) / mj * R_quad * e_3 - Xdd_des_quad;
        %e_pdd_hat = (Xdd - Xdd_des(:, i))/dt_sim;
        
        nu_ej = force_prev(j) * R_quad * e_3 /mj - Xdd_des_quad + kv_j * e_pd / mj ...
                + kp_j * e_p / mj + gravity + delta_hat(:, j) / mj;
        eta = - alpha * mj * nu_ej ...
              - mj * gamma * (e_pd + epsilon * e_p);
        vec = R_quad' * (eta - mj * Xddd_des_quad(j) + kv_j * e_pdd_hat + kp_j * e_pd + delta_hat(:, j));
        w_xj_des = - vec(2) / force_prev(j);
        w_yj_des = vec(1) / force_prev(j);
        forced = vec(3);

        force_j = force_prev(j) + forced * dt_sim * delay_bs;
        force_prev(j) = force_j;

        w_quad_des = [w_xj_des; w_yj_des; 0.0]; 
        end
               
        % rotaion control
        if mod(i, delay_quad) == 1 || delay_quad == 1
        Ij = It + Ib;
        kw_j = kw_I * norm(Ij);
        
        R_e_j = R * R_shape{j};
        w_e_j = R_shape{j}' * w;


        R_e_des_j = R_e_des{i} * R_shape{j};
        w_e_des_j = R_shape{j}' * w_e_des(:, i);
        
        % End effector control
        pitch_des = atan2(R_e_des_j(1,3), R_e_des_j(1,1));
        pitch = atan2(R_e_j(1,3), R_e_j(1,1));
        e_pitch = wrapToPi(pitch - pitch_des);
        w_e_des_j = w_e_des_j(2);
        w_e_ref = w_e_des_j - k_pitch * e_pitch;
        thetad_ref = w_e_ref - w_quad_des(2);
        
        if i > 1
            wd_quad_des = (w_quad_des - w_quad_des_prev(:, j)) / dt_sim / delay_quad;
        else
            wd_quad_des = [0; 0; 0];
        end
        w_quad_des_prev(:, j) = w_quad_des;
        
        e_w_quad = w_quad_j - w_quad_des;
        
        %quad rotor : kinematic level?
        tau(j) = (It(2,2) + Ib(2,2)) * wd_quad_des(2)  - kw_j * e_w_quad(2); 
        %tau(j) = tau(j) - disturb(2)/num_AMs;
        end

        % servo motor
        theta_ref(j) = theta_ref(j) + thetad_ref * dt_sim;
        e_theta = theta(j) - theta_ref(j) ;
        e_thetad = thetad(j) - thetad_ref;
        tau_theta(j) = - kp_servo * e_theta - kd_servo * e_thetad - damp_servo * thetad(j); % 1-dim

        % TODO
        %tau(j) = tau(j) + tau_theta(j);
        %tau(j) = tau(j) - disturb(2)/num_AMs;
 
        thrust_j = inv(B) * [0; tau(j); 0; force_j];
        thrust_j = min(thrust_limit, max(-thrust_limit, thrust_j));
        
        tau(j) = B(2,:) * thrust_j;
        force_j = B(4, :) * thrust_j;
        force_j = force_j; 

        term = S(rj) * R_quad * [0; 0; force_j];
        tau_tot(2) = tau_tot(2) + tau_theta(j) + term(2); % 1-dim
        force_tot = force_tot + force_j * Rt{j} * e_3; % 3 - dim
        
        phidd = It(2,2) \ (tau(j) - tau_theta(j));
        phid(j) = phid(j) + phidd * dt_sim;
        phi(j) = phi(j) + phid(j) * dt_sim + 0.5 * phidd * dt_sim^2 ;
        Rt{j} = Ry(phi(j));

        theta(j) = wrapToPi(atan2(R_e_j(1,3), R_e_j(1,1)) - phi(j));
        thetad_prev = thetad(j);
        thetad(j) = wrapToPi(w_e_j(2) - phid(j));
        thetadd(j) = (thetad(j) - thetad_prev) / dt_sim;

        thrusts = [thrusts; thrust_j];
        force_per_M = [force_per_M; force_j / mj];
        delta_hat_x_per_M = [delta_hat_x_per_M; delta_hat(3, j) / mj];
    end
    % generate disturbance
    disturb_dot = randn(6, 1) * sigma;
    disturb = disturb + disturb_dot * dt_sim;
    disturb = min(max(disturb, - max_val), max_val);
    disturb_sim(i,:) = disturb;

    %force_tot : sigma lambdai *Ri *e3
    Xdd = AM_mass\(force_tot - AM_mass * 9.81 * e_3 + disturb(4:6));
    % tau_tot : sigma tau_i + tau_theta(1) + r x lambda
    wd = inv(AM_inertia) * (tau_tot - S(w) * AM_inertia * w + [0;disturb(2);0]);

    Xd = Xd + Xdd * dt_sim;
    X = X + Xd * dt_sim + 0.5 * Xdd * dt_sim^2 ;
    w = w + wd * dt_sim;
    Rd = R * S(w);
    R = R + Rd * dt_sim;
    [U_, ~, V_] = svd(R);
    R  = U_ * V_';

    % Log
    X_hist  = [X_hist, X];
    Xd_hist = [Xd_hist, Xd];
    w_hist  = [w_hist, w_quad_j];
    w_des_hist = [w_des_hist, w_quad_des];
    wd_hist = [wd_hist, wd];
    wd_des_hist = [wd_des_hist, wd_quad_des]; 
    R_hist{i} = R;
    thrusts_hist = [thrusts_hist, thrusts];
    force_per_M_hist = [force_per_M_hist, force_per_M];
    delta_hat_x_per_M_hist = [delta_hat_x_per_M_hist, delta_hat_x_per_M];
    
    % totoal
    e_p = X - X_des(:, i); 
    e_pd = Xd - Xd_des(:, i);
    
    % quad
    e_w_quad = [0; phid(j); 0] - w_quad_des;
   
    e_p_hist = [e_p_hist, e_p]; 
    e_d_hist = [e_d_hist, e_pd]; 
    e_w_hist = [e_w_hist, e_w_quad];
    nu_e_hist = [nu_e_hist, nu_ej];
    delta_hat_hist = [delta_hat_hist, delta_hat(:, j)];
    delta_tilde_hist = [delta_tilde_hist, mj / AM_mass * disturb(4:6) - delta_hat(:, j)];
    
    pitch_hist = [pitch_hist, wrapToPi(pitch)];
    pitch_des_hist = [pitch_des_hist, wrapToPi(pitch_des)];
    phi_hist = [phi_hist, phi];
   
    e_pitch_hist = [e_pitch_hist, e_pitch];
    e_theta_hist = [e_theta_hist, e_theta];
    e_thetad_hist = [e_thetad_hist, e_thetad]; 

    times = [times; i * dt_sim];
end
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
    plot(times, w_des_hist(j, 1:i), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
end
legend({'$\omega_y$', '$\omega_y^{\mathrm{des}}$'}, ...
       'Interpreter', 'latex','FontSize', 12);
title('${\omega}$ vs Desired', 'Interpreter', 'latex','FontSize', 14)
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
    plot(times, delta_hat_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$\hat{\Delta}_{p,i}$', 'Interpreter', 'latex','FontSize', 14)
grid on

%delta_tilde
subplot(3,2,6)
hold on
for j = 1:3
    plot(times, delta_tilde_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
end
legend({'$x$','$y$','$z$'},'Interpreter','latex','FontSize', 12);
title('$\tilde{\Delta}_{p,i}$', 'Interpreter', 'latex','FontSize', 14)
grid on
%% Assumption: force/disturbance balance
figure('Position',[1350 550 500 400]);

% 7. thrusts
subplot(3,1,1)
hold on
plot(times, thrusts_hist, 'LineWidth', 1.0);
%ylim([ -thrust_limit, thrust_limit] )
title('Thrusts', 'Interpreter', 'latex','FontSize', 14)
grid on


legend_entries = arrayfun(@(x) sprintf('%d', x), 1:num_AMs, 'UniformOutput', false);
subplot(3,1,2)
plot(times, force_per_M_hist)
title('$\frac{\lambda_{p,i}}{M_i}$', 'Interpreter', 'latex','FontSize', 14)
legend(legend_entries)
grid on

subplot(3,1,3)
plot(times, delta_hat_x_per_M_hist)
title('$\frac{\Delta_{p,i}}{M_i}$', 'Interpreter', 'latex','FontSize', 14)
legend(legend_entries)
grid on
%% End effector Error 
figure('Position',[1350 50 500 450]);
colors = lines(6);

subplot(3,1,1)
hold on
for j = 1:1
    plot(times, e_pitch_hist(j, :) / pi * 180, 'Color', colors(3-j,:), 'LineWidth', 1.0);
end
ylabel("degree")
legend({'$pitch$'},'Interpreter','latex','FontSize', 12);
title('$e_{pitch}$', 'Interpreter', 'latex','FontSize', 14)
grid on


subplot(3,1,2)
hold on
for j = 1:1
    plot(times, e_theta_hist(j, :), 'Color', colors(3-j,:), 'LineWidth', 1.0);
end
legend({'$pitch$'},'Interpreter','latex','FontSize', 12);
title('$e_\theta$', 'Interpreter', 'latex','FontSize', 14)
grid on

subplot(3,1,3)
hold on
for j = 1:1
    plot(times, e_thetad_hist(j, :), 'Color', colors(3-j,:), 'LineWidth', 1.0);
end
legend({'$pitch$'},'Interpreter','latex','FontSize', 12);
title('$e_{\dot\theta}$', 'Interpreter', 'latex','FontSize', 14)
grid on
%% Video
figure('Position',[600 100 800 800]);
dN = 0.2 / dt_sim;
framesPerSecond = 1/dt_sim/dN;
rate = rateControl(framesPerSecond);
arrow_len = 0.2;
for i = 1:dN:N_sim_tmp
    clf;
    grid on; axis equal;
    xlim([X_des(1, i) - (l1+l2)/2*num_AMs*1.3, X_des(1, i) + (l1+l2)/2*num_AMs*1.3 ] )
    ylim([X_des(3, i) - (l1+l2)/2*num_AMs*1.3, X_des(3, i) + (l1+l2)/2*num_AMs*1.3 ])
    
    hold on;
    pitch = pitch_hist(i);
    R = Ry(pitch);
    
    for j = 1:num_AMs
        R_e_j = R * R_shape{j};
        X_quad = X_hist(:, i) + R * r_cj{j};
        R_quad = Ry(phi_hist(j, i));
        dx = force_per_M_hist(j, i) / 9.8 * arrow_len * R_quad(1,3);
        dz = force_per_M_hist(j, i) / 9.8 * arrow_len * R_quad(3,3);       
        plot(X_quad(1), X_quad(3), 'o', 'Color', 'b');
        quiver(X_quad(1), X_quad(3), dx, dz, 0, 'r', 'LineWidth', 1.5);
        
        X1 = X_quad + l1 * R_e_j * e_1;
        X2 = X_quad - l2 * R_e_j * e_1;
        plot([X1(1), X2(1)], [X1(3), X2(3)], 'b-', 'LineWidth', 1.0);
   
    end

    pitch_des = pitch_des_hist(i);
    R_des = Ry(pitch_des);
    
    for j = 1:num_AMs
        R_e_des_j = R_des * R_shape{j};
        X_des_quad = X_des(:, i) + R_des * r_cj{j};
        %dx = arrow_len * R_e_des_j(1,1);
        %dz = arrow_len * R_e_des_j(3,1);       
        plot(X_des_quad(1), X_des_quad(3), 'o', 'Color', 'black', 'LineWidth', 2.0);
        %quiver(X_des_quad(1), X_des_quad(3), dx, dz, 0, 'black', 'LineWidth', 2.0);

        X1 = X_des_quad + l1 * R_e_des_j * e_1;
        X2 = X_des_quad - l2 * R_e_des_j * e_1;
        plot([X1(1), X2(1)], [X1(3), X2(3)], 'k--', 'LineWidth', 2.0);
    end

    % com
    plot(X_hist(1, i), X_hist(3, i), 'o', 'Color', 'b', 'MarkerSize', 10);
    plot(X_des(1, i), X_des(3, i), 'o', 'Color', 'k', 'MarkerSize', 8);

    xlabel('X position'); ylabel('Z position');
    title_string = sprintf("time : %.2f sec", times(i));
    title(title_string);
    drawnow
    %waitfor(rate);
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