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
l1 = 0.35; l2 = 0.30; % ca = 0.24

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

M = zeros(num_AMs * 6, num_AMs * 6);
for j =1:num_AMs
    M(3*j-2: 3*j, 3*j-2: 3*j) = mass_ams(j) * eye(3);
    M(3*(j+num_AMs)-2: 3*(j+num_AMs), 3*(j+num_AMs)-2: 3*(j+num_AMs)) = It + Ib;
end
%%
wn = 2; damp = 1.5; % 2, 1.1
kp_M = wn^2; 
kv_M = 2 * damp *sqrt(kp_M);

wn = 2; damp = 1.1; % 2, 1.1
kp_z = wn^2; 
kv_z = 2 * damp *sqrt(kp_z);

kp_M = diag([kp_M, kp_M, kp_z]);
kv_M = diag([kv_M, kv_M, kv_z]);

kw_I = 10; % 10 or 20 or 100

% Backstepping gain
epsilon = kv_M(1, 1) * 0.3; % kv_M(1, 1) * 0.7
alpha = 4; % 4 or 10
gamma = alpha * 1.0; % alpha * 1.0 

% servo moter
kp_servo = 0.01; % 0.1 / 0.01
kd_servo = 0.1; % 1.0 / 0.1 / 0.03
damp_servo = 0.05; % 0.05

% endeffector control
damped = 0.0; % 0.3
k_pitch = 0.5; % 1.0

% MBO gain : 게인값을 높이고 필터링도 괜찮을듯 
k_mbo = diag([1.0 1.0 1.0 1.0]) * 3e0; % 3e0
dt_lpf = 0.2;

% optimization
tan_max = tan(45 / 180 * pi);
% k_internal_f = 1e0; % 1e0 ~, might be depend on external wrench
% k_internal_tau = 3e0;% 1e1 ~
k_smooth1 = 1e0;% 1e-2 ~ 1e0
dt_lpf_delta = 0.01;

f_int_max = 100;
tau_int_max = 100;


% Simulation Parmeters
N_sim_tmp = 5000;
show_video = true;
save_video = true;
video_speed = 1.0;

% Thrust limit and w_des limit
thrust_limit = thrust_limit * 1.0; % 1 ~ 3
w_des_limit = 3.0; % 2 ~ 1

% disturbance
% payload at EE
payload = 5;
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
delay_mbo = 0.01 / dt_sim; % 0.01

rng('shuffle')

X_hover = [1; 0; 3] * 1e-1;
rpy_hover = [0, 10, 0] / 180 * pi; velocity = [-0.0; 0; 0]; maximum_X = [10.0; 0; 10.0];
[X_des, Xd_des, Xdd_des, Xddd_des, R_e_des, w_e_des, wd_e_des] = get_traj_hover_manip(X_hover, rpy_hover, velocity, maximum_X, N_sim, dt_sim);
% helix
radius = 0.3;  v_z = 0.0;
omega = 2 * pi * 0.3; 
rpyd  = [0.00; -0.00; 0.0] * 2 * pi;
X_hover = [0.1; 0; 0.3]; rpy_hover = [0, 0, 0] / 180 * pi; 
%[X_des, Xd_des, Xdd_des, Xddd_des, R_e_des, w_e_des, wd_e_des] = get_traj_helix_manip_2d(radius, omega, v_z, rpyd, X_hover, rpy_hover, N_sim, dt_sim);
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

lambda_prev = mass_ams * norm(gravity) * 1.0;
delta_hat = zeros(4, num_AMs);
delta_hatd = zeros(4, num_AMs);
f_ext_hat = zeros(4, 1);
h_0 = [AM_mass * Xd; e_2' * (AM_inertia - It * num_AMs) * w + It(2,2) * sum(phid)]; % linear moment
int_mbo = zeros(4, 1);
force_tot_mbo = AM_mass * norm(gravity) * e_3;
tau_tot_mbo = zeros(3, 1);
lambda_mbo = lambda_prev;
phi_prev = phi; phid_prev = phid;
delta_prev = zeros(3*num_AMs, 1); delta_pprev = zeros(3*num_AMs, 1);

w_quad_des_prev = zeros(3, num_AMs);

X_hist = []; Xd_hist = []; w_hist = []; wd_hist = []; R_hist = cell(N_sim, 1); 
w_des_hist = []; wd_des_hist = []; thrusts_hist = []; F_hat_hist = [];
e_p_hist = []; e_d_hist = []; e_w_hist = []; e_psi_hist = [];
nu_e_hist = []; delta_hat_hist = []; delta_tilde_hist = [];
force_per_M_hist = []; delta_hat_x_per_M_hist = []; delta_hat_z_per_M_hist = [];
e_theta_hist = []; e_thetad_hist = [];
e_pitch_hist = []; e_R_hist = [];
pitch_hist = []; pitch_des_hist = []; phi_hist = [];
disturb_hist = []; f_ext_hat_hist = [];
internal_f_opt_hist = []; internal_tau_opt_hist = [];
internal_f_real_hist = []; internal_tau_real_hist = []; 
internal_f_real_tq_hist = []; internal_tau_real_tq_hist = []; 
internal_f_real_qd_hist = []; internal_tau_real_qd_hist = []; 
times = [];
R_hist = cell(N_sim_tmp, 1);

tic
for i = 1:N_sim_tmp
    if mod(i, N_sim_tmp/100) == 0
        fprintf("\n%d / %d\n", i, N_sim_tmp)
    end
    tau_tot = zeros(3, 1);
    if mod(i, delay_quad) == 1 || delay_quad == 1
        tau = zeros(1, num_AMs);
    end
    tau_theta = zeros(num_AMs, 1);
    force_tot = zeros(3, 1);
    thrusts = [];
    force_per_M = [];
    delta_hat_x_per_M = [];
    delta_hat_z_per_M = [];
    
    input_real = zeros(num_AMs*6, 1);

    % estimation error
    X_error_dot = randn(3, num_AMs) * sigma_X;
    X_error = X_error + X_error_dot * dt_sim;
    X_error = min(max(X_error, - max_X), max_X);
    w_error_dot = randn(3, num_AMs) * sigma_w;
    w_error = w_error + w_error_dot * dt_sim;
    w_error = min(max(w_error, - max_w), max_w);
    
    % MBO
    if mod(i, delay_mbo) == 1 || delay_mbo == 1
        h = [AM_mass * Xd; e_2' * (AM_inertia - It * num_AMs) * w + It(2,2) * sum(phid)];
        int_mbo = int_mbo + ([(force_tot_mbo - AM_mass * norm(gravity) * e_3); tau_tot_mbo(2)] + f_ext_hat) * dt_sim * delay_mbo; 
        % 1-order LPF
        f_ext_hat_sensored = k_mbo * (h - int_mbo - h_0);
        f_ext_hat = (dt_lpf * f_ext_hat + dt_sim * delay_mbo * f_ext_hat_sensored) / (dt_lpf + dt_sim * delay_mbo);
        
        tau_tot_mbo = zeros(3, 1);
        

        % Internal forces : expensive
        A = zeros((num_AMs-1)*6,num_AMs*6);          
        Ad = zeros((num_AMs-1)*6,num_AMs*6);          
        % l1_hat = l1 * kinematic_error;
        % l2_hat = l2 * kinematic_error;
        
        for j = 1:num_AMs - 1
            A(3*j-2: 3*j, 3*j-2: 3*j+3) = [eye(3), -eye(3)];
            A(3*j-2: 3*j, 3*(j+num_AMs)-2: 3*(j+num_AMs)+3) = [l2 * R * R_shape{j} * S(e_1), l1 * R * R_shape{j+1} * S(e_1)];
            A(3*(j+num_AMs-1)-2: 3*(j+num_AMs-1), 3*(j+num_AMs)-2: 3*(j+num_AMs)+3) = [eye(3), -eye(3)];
            Ad(3*j-2: 3*j, 3*(j+num_AMs)-2: 3*(j+num_AMs)+3) = [l2 * R * R_shape{j}  * S(w) * S(e_1), l1 * R * R_shape{j+1}  * S(w) * S(e_1)];
        end
        A_dagger = ((A* (M\ A')) \ A) * inv(M);
        
        % SOCP
        f = [zeros(3*num_AMs, 1); 1; 1];
        
        % External wrench compensation
        A_eq = zeros(6, 3*num_AMs + 2);
        r_e = R * r_cj{1} + l1 * R * R_shape{1} * e_1;
        for j = 1:num_AMs
            rj = R * r_cj{j} - r_e;
            A_eq(1:6, 3*j-2:3*j) = [eye(3,3); S(rj)];
        end
        b_eq = [f_ext_hat(1:3); e_2 * f_ext_hat(4) - S(r_e) * f_ext_hat(1:3)]; % w.r.t ee 
        
        % Each module's decentralized input
        v = zeros(3*num_AMs, 1);
        feedback = zeros(3*num_AMs, 1);
        for j = 1:num_AMs
            mj = mass_ams(j);
            rj = r_cj{j};
            kp_j = mj * kp_M;
            kv_j = mj * kv_M;
            
            X_quad = X + R * rj; 
            X_des_quad = X_des(:, i) + R_e_des{i} * rj;
            Xd_quad = Xd + Rd * rj; 
            Xd_des_quad = Xd_des(:, i) + R_e_des{i} * S(w_e_des(:, i)) * rj;
            Xdd_des_quad = Xdd_des(:, i) + R_e_des{i}* S2(w_e_des(:, i)) * rj + R_e_des{i}* S(wd_e_des(:,i)) * rj;

            e_p  = (1 - X_error(:, j)).*X_quad - X_des_quad; 
            e_pd = (1 - X_error(:, j)).*Xd_quad - Xd_des_quad;
            
            feedback(3*j-2:3*j) = mj * Xdd_des_quad - kv_j *e_pd - kp_j * e_p;
            v(3*j-2:3*j) = feedback(3*j-2:3*j) + mj * 9.81 * e_3; %TODO feed-back term
        end
        
        % inf-norm of v - delta
        socConstraints = [];
        for j = 1:num_AMs
            A_soc = [zeros(3, 3*(j-1)), eye(3), zeros(3, 3*(num_AMs-j)+2)];
            b_soc = [v(3*j-2:3*j)];
            d_soc = [zeros(3*num_AMs,1); 1; 0];
            gamma = 0;
            socConstraints = [socConstraints, secondordercone(A_soc, b_soc, d_soc, gamma)];
        end

        % tilting angle constraint
        A1 = zeros(num_AMs, 3*num_AMs);
        A2 = zeros(num_AMs, 3*num_AMs);
        a1 = - e_1' + tan_max * e_3';
        a2 = e_1' + tan_max * e_3';
        for j = 1 : num_AMs
            Rj = R * R_shape{j};
            A1(j, 3*j-2:3*j) = a1 * Rj';
            A2(j, 3*j-2:3*j) = a2 * Rj';
        end
        
        A_ineq = [A1, zeros(num_AMs, 2);...
                  A2, zeros(num_AMs, 2)];
        b_ineq = [A1 * v; A2 * v];
        
        % inf-norm of internal force % TODO: torque check
        r_e = R * r_cj{1} + l1 * R * R_shape{1} * e_1;
        F_ext = [f_ext_hat(1:3); zeros(3*num_AMs -3, 1);
                 e_2 * f_ext_hat(4) - S(r_e) * f_ext_hat(1:3); zeros(3*num_AMs -3, 1)];
        qd = [zeros(3*num_AMs, 1); repmat(w, num_AMs, 1)];
        Cqd = zeros(6*num_AMs) * qd;
        F_int0 = A_dagger *([v; zeros(3*num_AMs, 1)] - F_ext + Cqd) - ((A* (M\ A')) \ Ad) * qd; % mge_3 doesn't generate internal forces;
        A_delta = A_dagger(:, 1:3*num_AMs);
        F_int_max = [f_int_max * ones(3*(num_AMs-1),1); tau_int_max * ones(3*(num_AMs-1),1)];

        A_ineq = [A_ineq; A_delta, zeros(6*(num_AMs-1), 2);...
                          -A_delta, zeros(6*(num_AMs-1), 2)];
        b_ineq = [b_ineq; -F_int0 + F_int_max; F_int0 + F_int_max];
        
        % inf-norm of second derivative
        K_smooth1_inv = k_smooth1 \ ones(num_AMs*3, 1);
        
        A_ineq = [A_ineq; eye(3*num_AMs), zeros(3*num_AMs, 1), -K_smooth1_inv;...
                          -eye(3*num_AMs), zeros(3*num_AMs, 1), -K_smooth1_inv];
        ref = delta_prev;
        b_ineq = [b_ineq; ref; -ref];
        
        lb = -inf * ones(3* num_AMs + 2, 1);   
        ub = inf * ones(3* num_AMs + 2, 1);            

        options = optimoptions('coneprog', ...
            'Display', 'off', ...             % 반복 과정 표시
            'LinearSolver', 'auto');       % 희소 문제에 적합한 솔버 지정
        
        [sol, fval, exitflag, solver_output] = coneprog(f, socConstraints, A_ineq, b_ineq, A_eq, b_eq, lb, ub, options);
        
        if exitflag == 1
            delta_prev = sol(1:3*num_AMs);
            delta_pprev = delta_prev;
            delta_reshape = reshape(sol(1:3*num_AMs), 3, num_AMs);
            norms = sol(3*num_AMs+1 : end);
        else
            disp(solver_output)
        end
        
        delta_hatd_computed = ([delta_reshape; zeros(1, num_AMs)] - delta_hat) / dt_sim / delay_mbo;
        delta_hatd = (dt_lpf_delta * delta_hatd + dt_sim * delay_mbo *delta_hatd_computed) / (dt_lpf + dt_sim * delay_mbo);
        delta_hat = [delta_reshape; zeros(1, num_AMs)];

        internal_wrench = A_delta *  delta_prev + F_int0;
        internal_f_opt = internal_wrench(1:3*num_AMs-3);
        internal_tau_opt = internal_wrench(3*num_AMs-2 :end);
        
    end

    % for j = num_AMs:-1:1
    for j = 1:num_AMs
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
            
            %TODO: sign of delta
            e_pdd_hat = -9.81 * e_3 + (lambda_prev(j) * R_quad * e_3  + delta_hat(1:3, j)) / mj - Xdd_des_quad;
            %e_pdd_hat = (Xdd - Xdd_des(:, i))/dt_sim/delay_bs;
            
            nu_ej    = lambda_prev(j) * R_quad * e_3 /mj - Xdd_des_quad + kv_j * e_pd / mj ...
                    + kp_j * e_p / mj - 9.81 * e_3 + delta_hat(1:3, j) / mj;
            eta = - alpha * mj * nu_ej ...
                  - mj * gamma * (e_pd + epsilon * e_p);
            vec = R_quad' * (eta + mj * Xddd_des_quad(:, j) - kv_j * e_pdd_hat - kp_j * e_pd - delta_hatd(1:3, j)); % delta_hatd
            w_xj_des = - vec(2) / lambda_prev(j);
            w_yj_des = vec(1) / lambda_prev(j);
            lambdad = vec(3);
    
            lambda_prev(j) = lambda_prev(j) + lambdad * dt_sim * delay_bs;
            
            w_quad_des = [w_xj_des; w_yj_des; 0.0]; 
            w_quad_des = max(-w_des_limit, min(w_des_limit, w_quad_des)); 
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
            
            % quad rotor : kinematic level?
            tau(j) = (It(2,2) + Ib(2,2)) * wd_quad_des(2) - kw_j * e_w_quad(2);
            tau(j) = tau(j) - delta_hat(4, j);
        end

        % servo motor
        theta_ref(j) = theta_ref(j) + thetad_ref * dt_sim;
        e_theta = theta(j) - theta_ref(j) ;
        e_thetad = thetad(j) - thetad_ref;
        tau_theta(j) = - kp_servo * e_theta - kd_servo * e_thetad - damp_servo * thetad(j); % 1-dim

        % saturated thrusts
        thrust_j = inv(B) * [0; tau(j); 0; lambda_prev(j)];
        thrust_j = min(thrust_limit, max(-thrust_limit, thrust_j));
        
        tau_j = B(2,:) * thrust_j;
        lambda_j = B(4, :) * thrust_j;

        term = S(rj) * R' * R_quad * [0; 0; lambda_j];
        tau_tot(2) = tau_tot(2) + tau_theta(j) + term(2); % 1-dim
        force_tot = force_tot + lambda_j * R_quad * e_3; % 3 - dim

        if mod(i, delay_mbo) == 1 || delay_mbo == 1
            tau_tot_mbo(2) = tau_tot_mbo(2) + tau_j + term(2);
        end
        
        phidd = It(2,2) \ (tau_j - tau_theta(j));
        phid(j) = phid(j) + phidd * dt_sim;
        phi(j) = phi(j) + phid(j) * dt_sim + 0.5 * phidd * dt_sim^2 ;
        Rt{j} = Ry(phi(j));

        theta(j) = wrapToPi(atan2(R_e_j(1,3), R_e_j(1,1)) - phi(j));
        thetad_prev = thetad(j);
        thetad(j) = wrapToPi(w_e_j(2) - phid(j));
        thetadd(j) = (thetad(j) - thetad_prev) / dt_sim;
        
        input_real(3*j -2:3*j) = R_quad * [0; 0; lambda_j] + mass_ams(j) * gravity;
        input_real(3*(j+num_AMs) -2:3*(j+num_AMs)) = (tau_j - tau_theta(j)) * e_2;
        thrusts = [thrusts; thrust_j];
        force_per_M = [force_per_M; lambda_j / mj];
        delta_hat_x_per_M = [delta_hat_x_per_M; delta_hat(1, j) / mj];
        delta_hat_z_per_M = [delta_hat_z_per_M; delta_hat(3, j) / mj];
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

    %force_tot : sigma lambdai *Ri *e3
    Xdd = (AM_mass * mass_uncertainty)\(force_tot - (AM_mass * mass_uncertainty) * 9.81 * e_3 + disturb(1:3));
    % tau_tot : sigma tau_i + tau_theta(1) + r x lambda
    wd = inv(AM_inertia * inertia_uncertainty - It * num_AMs) * (tau_tot - S(w) * (AM_inertia * inertia_uncertainty - It * num_AMs) * w + disturb(4:6));
    
    if mod(i, delay_mbo) == 1 || delay_mbo == 1
        force_tot_mbo = force_tot;
        phi_prev = phi;
        phid_prev = phid;
        lambda_mbo = lambda_prev;
    end

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
    delta_hat_z_per_M_hist = [delta_hat_z_per_M_hist, delta_hat_z_per_M];
    
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
    delta_tilde_hist = [delta_tilde_hist, mj / AM_mass * disturb(1:3) - delta_hat(1:3, j)];
    f_ext_hat_hist = [f_ext_hat_hist, f_ext_hat];
    disturb_hist = [disturb_hist, disturb];

    pitch_hist = [pitch_hist, wrapToPi(pitch)];
    pitch_des_hist = [pitch_des_hist, wrapToPi(pitch_des)];
    phi_hist = [phi_hist, phi];
   
    e_pitch_hist = [e_pitch_hist, e_pitch];
    e_theta_hist = [e_theta_hist, e_theta];
    e_thetad_hist = [e_thetad_hist, e_thetad]; 
    times = [times; i * dt_sim];
    
    % Internal force log
    internal_f_opt_hist = [internal_f_opt_hist, internal_f_opt]; 
    internal_tau_opt_hist = [internal_tau_opt_hist, internal_tau_opt];

    F_ext = [f_ext_hat(1:3); zeros(3*num_AMs -3, 1);
    e_2 * f_ext_hat(4) - S(r_e) * f_ext_hat(1:3); zeros(3*num_AMs -3, 1)];
    qd = [zeros(3*num_AMs, 1); repmat(w, num_AMs, 1)];
    Cqd = zeros(6*num_AMs) * qd;

    internal_qd = - ((A* (M\ A')) \ Ad) * qd;
    internal_f_real_qd_hist = [internal_f_real_qd_hist, internal_qd(1:3*(num_AMs-1))]; 
    internal_tau_real_qd_hist = [internal_tau_real_qd_hist, internal_qd(3*(num_AMs-1)+1:end)];
    
    internal_tq = A_dagger * [zeros(3*num_AMs, 1); input_real(3*num_AMs+1:end)];
    internal_f_real_tq_hist = [internal_f_real_tq_hist, internal_tq(1:3*(num_AMs-1))]; 
    internal_tau_real_tq_hist = [internal_tau_real_tq_hist,  internal_tq(3*(num_AMs-1)+1:end)];

    internal_real = A_dagger * (-input_real - F_ext - Cqd) + internal_qd;
    internal_f_real_hist = [internal_f_real_hist, internal_real(1:3*(num_AMs-1))]; 
    internal_tau_real_hist = [internal_tau_real_hist, internal_real(3*(num_AMs-1) + 1:end)];
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
    plot(times, w_des_hist(j, 1:i), '--', 'Color', colors(j,:), 'LineWidth', 2.5);
end
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
for j = 1:4
    plot(times, f_ext_hat_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0);
    if j==4 
        plot(times, disturb_hist(5, :), 'Color', colors(j,:), 'LineWidth', 1.0, 'LineStyle','--');
    else
        plot(times, disturb_hist(j, :), 'Color', colors(j,:), 'LineWidth', 1.0, 'LineStyle','--');
    end
end
legend({'$\hat{\Delta}_x^p$','$\Delta_x^p$',...
        '$\hat{\Delta}_y^p$','$\Delta_y^p$',...
        '$\hat{\Delta}_z^p$','$\Delta_z^p$',...
        '$\hat{\Delta}_y^\omega$','$\Delta_y^\omega$'},'Interpreter','latex','FontSize', 12);
title('$\hat{\Delta}$', 'Interpreter', 'latex','FontSize', 14)
grid on

%delta_tilde
subplot(3,2,6)
hold on
for j = 1:4
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

%% force/disturbance
figure('Position',[1350 300 400 600]);

% 7. thrusts
subplot(4,1,1)
hold on
plot(times, thrusts_hist, 'LineWidth', 1.0);
ylabel("[N]")
ylim([ 0, thrust_limit*1] )
title('Thrusts', 'Interpreter', 'latex','FontSize', 14)
grid on


legend_entries = arrayfun(@(x) sprintf('%d', x), 1:num_AMs, 'UniformOutput', false);
subplot(4,1,2)
plot(times, force_per_M_hist * mass_ams(1))
title('$\lambda_{i}$', 'Interpreter', 'latex','FontSize', 14)
ylabel("[N]")
ylim([ 0, thrust_limit*4] )
legend(legend_entries)
grid on

subplot(4,1,3)
plot(times, delta_hat_x_per_M_hist)
title('$\frac{\Delta_{p,i}}{M_i} x $', 'Interpreter', 'latex','FontSize', 14)
legend(legend_entries)
grid on

subplot(4,1,4)
plot(times, delta_hat_z_per_M_hist)
title('$\frac{\Delta_{p,i}}{M_i} z $', 'Interpreter', 'latex','FontSize', 14)
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
%% Internal forces
figure('Position',[1000 50 600 700]);
colors = lines(num_AMs*3);

subplot(4,2,1)
hold on
for j = 1:num_AMs-1
    % plot(times, internal_f_opt_hist(3*j-2:3*j, :), 'LineWidth', 1.0);
    plot(times, internal_f_opt_hist(3*j-2, :), 'LineWidth', 1.0);
    plot(times, internal_f_opt_hist(3*j, :), 'LineWidth', 1.0);
end
% legend({'$2,x$','$2,z$','$3,x$','$2,z$' },'Interpreter','latex','FontSize', 8);
ylabel("[N]")
title('Internal forces(Optimized)', 'Interpreter', 'latex','FontSize', 14)
grid on


subplot(4,2,2)
hold on
for j = 1:num_AMs-1
    plot(times, internal_tau_opt_hist(3*j-1, :), 'LineWidth', 1.0);
end
ylabel("[Nm]")
title('Internal torques(Optimized)', 'Interpreter', 'latex','FontSize', 14)
grid on


subplot(4,2,3)
hold on
for j = 1:num_AMs-1
    % plot(times, internal_f_opt_hist(3*j-2:3*j, :), 'LineWidth', 1.0);
    plot(times, internal_f_real_hist(3*j-2, :), 'LineWidth', 1.0);
    plot(times, internal_f_real_hist(3*j, :), 'LineWidth', 1.0);
end
ylabel("[N]")
title('Internal forces(Real)', 'Interpreter', 'latex','FontSize', 14)
grid on

subplot(4,2,4)
hold on
for j = 1:num_AMs-1
    plot(times, internal_tau_real_hist(3*j-1, :), 'LineWidth', 1.0);
end
ylabel("[Nm]")
title('Internal torques(Real)', 'Interpreter', 'latex','FontSize', 14)
grid on

subplot(4,2,5)
hold on
for j = 1:num_AMs-1
    % plot(times, internal_f_opt_hist(3*j-2:3*j, :), 'LineWidth', 1.0);
    plot(times, internal_f_real_tq_hist(3*j-2, :), 'LineWidth', 1.0);
    plot(times, internal_f_real_tq_hist(3*j, :), 'LineWidth', 1.0);
end
ylabel("[N]")
title('Internal forces(by torque)', 'Interpreter', 'latex','FontSize', 14)
grid on

subplot(4,2,6)
hold on
for j = 1:num_AMs-1
    plot(times, internal_tau_real_tq_hist(3*j-1, :), 'LineWidth', 1.0);
end
ylabel("[Nm]")
title('Internal torques(by torque)', 'Interpreter', 'latex','FontSize', 14)
grid on

subplot(4,2,7)
hold on
for j = 1:num_AMs-1
    plot(times, internal_f_real_qd_hist(3*j-2:3*j, :), 'LineWidth', 1.0);
end
ylabel("[N]")
title('Internal forces(by Ad qd)', 'Interpreter', 'latex','FontSize', 14)
grid on

subplot(4,2,8)
hold on
for j = 1:num_AMs-1
    plot(times, internal_tau_real_qd_hist(3*j-1, :), 'LineWidth', 1.0);
end
ylabel("[Nm]")
title('Internal torques(by Ad qd)', 'Interpreter', 'latex','FontSize', 14)
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

for i = 1:dN:N_sim_tmp
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
        dx = force_per_M_hist(j, i) * mass_ams(j) / 9.81 * arrow_len * R_quad(1,3);
        dz = (force_per_M_hist(j, i) - 0 )* mass_ams(j) / 9.81 * arrow_len * R_quad(3,3);       
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