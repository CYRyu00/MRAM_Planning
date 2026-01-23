clear; close all;
addpath("functions", "../params",  "../../casadi-3.6.7-windows64-matlab2018b")
import casadi.*
%% Define Dynamic parameters
params = define_params_ver2();
mu = params{3}; r = params{4}; d = params{5};
thrust_limit= params{6}; gravity = params{16};
mb = params{17}; cb = params{18}; Ib = params{19};
mt = params{20}; ct = params{21}; It = params{22};
ma = params{23}; ca = params{24}; Ia = params{25};

% totoal robotic arm
Ib = Ia + Ib;% Ib = Ib * 0.01;
mb = ma + mb;
m0 = mt + mb;

B = [r r -r -r;  -r r r -r; mu -mu mu -mu; 1 1 1 1];

e_1 = [1; 0; 0]; e_2 = [0; 1; 0]; e_3 = [0; 0; 1];

mass_obj = 20;

box_width = 1.0; box_height = 0.6;
mu_s = 0.15;
mu_d = 0.10; 

%simulation params
v_s = 0.1; % should be v_s > mu_s * 9.81 *dt : a_max * dt
sigma = 0.0;

dt = 0.05; N = 10/dt;
n_am = 4;
mass_am = m0 * ones(1,n_am); m_com = sum(mass_am);
l1 = 0.35; l2 = 0.3;

x_i = [0; 0];
x_f = [-2.0; 0];
tan_max = tan(45 / 180 *pi);
f_int_max = 30;
tau_int_max = 3;
k_smooth = 1e1;
k_damp = 1e2;

% solver options
nlp_opts = struct;
nlp_opts.ipopt.print_level = 5;
nlp_opts.ipopt.tol = 1e-3;
nlp_opts.ipopt.max_iter = 1000;
nlp_opts.ipopt.mu_strategy = 'monotone'; % monotone, adaptive
nlp_opts.ipopt.linear_solver = 'mumps';
nlp_opts.ipopt.jacobian_approximation = 'exact'; % exact, limited-memory
nlp_opts.print_time = 5;

waypoints = [x_i(1), x_f(1)]; 
timePoints = [0, N*dt]; 
t_samples = (0:N)*dt;
velocityBoundaryCondition = [x_i(2), x_f(2)]; 
[xi_init, xid_init] = cubicpolytraj(waypoints, timePoints, t_samples, ...
                                 'VelocityBoundaryCondition', velocityBoundaryCondition);
x_n_init = [xi_init; xid_init];

% figure;
% hold on
% plot(t_samples, xi_init);
% plot(t_samples, xid_init);
% grid on
% title("x_n init")
% legend("xi", "xid")
%%
M_am = zeros(n_am * 6, n_am * 6);
for j =1:n_am
    M_am(3*j-2: 3*j, 3*j-2: 3*j) = mass_am(j) * eye(3);
    M_am(3*(j+n_am)-2: 3*(j+n_am), 3*(j+n_am)-2: 3*(j+n_am)) = It + Ib;
end
M_o = diag(mass_obj * ones(3,1));
M = [M_o, zeros(3, n_am*6); zeros(n_am*6, 3), M_am];

shape_pitch =[0 -5 -10 -15 -20]; % = [0 -10 -20 -20 -10 0 10] / 180 *pi;

AM_inertia = zeros(3, 3); % inertia w.r.t. its com
AM_com = [0; 0; 0];% 0 to com
r_0j = cell(n_am, 1); % 0 to j'th module
r_cj = cell(n_am, 1); % com to j'th module
I_cj = cell(n_am, 1); % inertia of j'th module w.r.t. com of shpe
R_shape = cell(1, n_am);
for j = 1:length(shape_pitch)
    R_shape{j} = Ry(shape_pitch(j) / 180 *pi);
end

AM_mass = sum(mass_am);
for j = 1:n_am
    if j == 1
        r_0j{j} = [0; 0; 0];
    else
        r_0j{j} = r_0j{j-1} -l2 * R_shape{j-1} * e_1 - l1 * R_shape{j} * e_1;% p_core to j
    end
    AM_com = AM_com + r_0j{j} * mass_am(j)/AM_mass;
end

% compute AM_inertia
for j = 1:n_am
    r_cj{j} = r_0j{j} - AM_com;% p_com to j
    I_cj{j} = Ib + It + mass_am(j) * (r_cj{j}' * r_cj{j} * eye(3,3) - r_cj{j} * r_cj{j}'); % TODO: compute Ib It for 3D 
    AM_inertia = AM_inertia + I_cj{j};
end

r_e_ = r_cj{1} + l1 * R_shape{1} * e_1; % position of EE w.r.t AM's com
r_o = r_e_ + [box_width/2; 0; - box_height/2]; % position of Object's com w.r.t AM's com
%% Dynamics
x = MX.sym('x', 2); % [ xi; xid ]
u = MX.sym('u', n_am * 3, 1); % gravity compensated force lamada R e3 = u + mge3, u = -delta

u_net = zeros(3, 1);
for j =1:n_am
    u_net = u_net + u(3*j-2:3*j);
end

%TODO
% 4th order runge-kutta
F_fric = -mass_obj * 9.81* (mu_d * tanh(4 * x(2) / v_s) ...
       + (mu_s - mu_d) *x(2) /v_s /((x(2)/2/v_s)^2 + 0.75)^2 ) - sigma * x(2);
k1 = [x(2); (e_1'* u_net + F_fric) / (mass_obj + m_com)];

x1 = x + k1 * dt/2;
F_fric = -mass_obj * 9.81* (mu_d * tanh(4 * x1(2) / v_s) ...
       + (mu_s - mu_d) *x1(2) /v_s /((x1(2)/2/v_s)^2 + 0.75)^2 ) - sigma * x(2);
k2 = [x1(2); (e_1'* u_net + F_fric) / (mass_obj + m_com)];

x2 = x + k2 * dt/2;
F_fric = -mass_obj * 9.81* (mu_d * tanh(4 * x2(2) / v_s) ...
       + (mu_s - mu_d) * x2(2) /v_s /((x2(2)/2/v_s)^2 + 0.75)^2 ) - sigma * x(2);
k3 = [x2(2); (e_1'* u_net + F_fric) / (mass_obj + m_com)];

x3 = x + k3 * dt;
F_fric = -mass_obj * 9.81* (mu_d * tanh(4 * x3(2) / v_s) ...
       + (mu_s - mu_d) * x3(2) /v_s /((x3(2)/2/v_s)^2 + 0.75)^2 ) - sigma * x(2);
k4 = [x3(2); (e_1'* u_net + F_fric) / (mass_obj + m_com)];

xd = (k1 + 2*k2 + 2*k3 + k4) /6;
x_next = x + xd * dt ;

% xd = [x(2); (e_1'* u_net - k_obj * x(1) - b_obj * x(2) )  / (mass_obj + m_com)];

% 1st order euler
% F_fric = -mass_obj * 9.81* (mu_d * tanh(4 * x(2) / v_s) ...
%        + (mu_s - mu_d) *x(2) /v_s /((x(2)/2/v_s)^2 + 0.75)^2 ) - sigma * x(2);
% xd = [x(2); (e_1'* u_net + F_fric) / (mass_obj + m_com)];
% x_next = x + xd * dt;

fd = Function('fd', {x, u}, {x_next});

% internal wrench 
A11_inv = zeros(n_am * 3, n_am * 3);
A12 = zeros(n_am * 3, n_am * 3);
A22_inv = zeros(n_am * 3, n_am * 3);
for j = 1:n_am
    % A11_inv(3*j - 2: 3*j, 3*j-2:end) = repmat(eye(3), 1, n_am - j +1);
    for k = j : n_am
        A11_inv(3*j - 2: 3*j, 3*k-2:3*k) = eye(3);
        A22_inv(3*j - 2: 3*j, 3*k-2:3*k) = Ry( (-shape_pitch(j) + shape_pitch(k)) /180*pi);
    end
    Rj = eye(3) * Ry(shape_pitch(j)/180*pi);
    if j < n_am
        A12(3*j - 2: 3*j, 3*j - 2: 3*j+3) = [ l1 * S(e_1) * Rj', l2 * S(e_1) * Rj' ];
    else
        A12(3*j - 2: 3*j, 3*j - 2: 3*j) = l1 * S(e_1) * Rj';
    end
end

qdd = MX.zeros(n_am*6, 1);
Cqd = MX.zeros(n_am*6, 1); % TODO for 3D
for j = 1:n_am
    qdd(6*j-5: 6*j) = [xd(2); 0; 0; 0; 0; 0]; % TODO pdd_com + Rcdd * r_ci
end
U = [u; zeros(n_am * 3, 1)];
AF_int = M_am * qdd + Cqd - U;
f_int_ = A11_inv * AF_int(1: n_am*3);
tau_int_ = - A22_inv * A12 * A11_inv * AF_int(1: n_am*3) + A22_inv * AF_int(n_am * 3 + 1: end);
f_int = Function('f_int', {x, u}, {f_int_});
tau_int = Function('tau_int', {x, u}, {tau_int_});

delta_n = MX.sym('delta_n', n_am * 3, N); % norminal delta = -u
x_n = MX.sym('x_n', 2, N + 1); % norminal trajectory [xi; pd_1]
t = MX.sym('t', 2, N ); % inf norm of 2-norm of thrusts, smoothing

opt_variables = [reshape(t, 2 * N, 1); reshape(x_n, 2 * (N + 1), 1); reshape(delta_n, n_am * 3 * N, 1)];
opt_var0 = [ones(1 * N, 1) * (mass_am(1) * 9.81)^2; zeros(1 * N, 1);
           reshape(x_n_init, 2 * (N + 1), 1); % x_n_init
            zeros(n_am * 3 * N, 1)]; % delta_n_init
%% Define NLP

% start position, end position
g = [x_n(:, 1); x_n(:, end)];
LBG = [x_i; x_f]; 
UBG = [x_i; x_f];

g = [g; delta_n(:, 1)];
LBG = [LBG; zeros(3*n_am, 1)];
UBG = [UBG; zeros(3*n_am, 1)];

obj = 0;
for i = 1: N
    %tiling angle ineq constraints
    for j = 1:n_am
        Rj = eye(3) * Ry(shape_pitch(j)/180*pi);
        feedforward = mass_am(j) * 9.81 * e_3 - delta_n(3*j-2:3*j, i);
        g = [g; (Rj*e_1)' * feedforward - tan_max * (Rj*e_3)' * feedforward;
               -(Rj*e_1)' * feedforward - tan_max * (Rj*e_3)' * feedforward];
        LBG = [LBG; -inf; -inf];
        UBG = [UBG; 0; 0];
    end
    %internal wrench constraints
    g = [g; f_int(x_n(:, i) , -delta_n(:, i)); ...
            tau_int(x_n(:, i) , -delta_n(:, i))];
    LBG = [LBG; -f_int_max * ones(3*n_am, 1); -tau_int_max * ones(3*n_am, 1)]; 
    UBG = [UBG; f_int_max * ones(3*n_am, 1); tau_int_max * ones(3*n_am, 1)];
    
    % Dynamics inequality
    g = [g; x_n(:, i+1) - fd(x_n(:, i) , -delta_n(:, i)) ];
    LBG = [LBG; zeros(2, 1)];
    UBG = [UBG; zeros(2, 1)];

    % object funtion
    if i > 1
        obj = obj + t(1, i);
        for j = 1:n_am
            feedforward = mass_am(j) * 9.81 * e_3 - delta_n(3*j-2:3*j, i);
            g = [g; feedforward' * feedforward - t(1, i)];
            LBG = [LBG; -inf];
            UBG = [UBG; 0];
        end
    end

    % if i > 1
    %     obj = obj + t(2, i);
    %     g = [g; k_smooth * (delta_n(:, i)- delta_n(:, i-1)) - t(2, i) * ones(3*n_am,1);
    %            -k_smooth * (delta_n(:, i)- delta_n(:, i-1)) - t(2, i) * ones(3*n_am,1)];
    %     LBG = [LBG; - inf * ones(6*n_am, 1)];
    %     UBG = [UBG; zeros(6*n_am, 1)];
    % end
    % 
    if i > 1
        obj = obj + k_smooth * (delta_n(:, i)- delta_n(:, i-1))' * (delta_n(:, i)- delta_n(:, i-1)) ;
    end

    obj = obj + k_damp * x_n(2, i)^2;


end

% solver options
nlp_prob = struct('x', opt_variables, 'f', obj, 'g', g);

solver = nlpsol('solver', 'ipopt', nlp_prob, nlp_opts);
sol = solver('x0', opt_var0, ...
             'lbx', -inf, 'ubx', inf,...
             'lbg', LBG, 'ubg', UBG);
fprintf("\nsuccess : %d\n", solver.stats.success)

solution = full(sol.x);
t_opt = reshape(solution(1:2*N), 2, N); % inf norm of 2-norm of thrusts, smoothing
x_n_opt = reshape(solution(2*N+1: 2*N + 2*(N + 1)), 2, N+1); % norminal trajectory [xi; xid]
delta_n_opt = reshape(solution(2*N + 2*(N + 1) + 1: end), 3*n_am, N); % norminal delta = -u

f_int_opt = zeros(3 * n_am, N);
tau_int_opt = zeros(3 * n_am, N);
for i = 1:N
    f_int_opt(:, i) = full(f_int(x_n_opt(:, i) , -delta_n_opt(:, i)));
    tau_int_opt(:, i) = full(tau_int(x_n_opt(:, i) , -delta_n_opt(:, i)));
end

%% trajectory save
% mkdir data
filename = sprintf('data/10sec_2m.mat');
N_opt = N; dt_opt = dt; n_am_opt = n_am; shape_pitch_opt = shape_pitch;
save(filename, 't_opt', 'x_n_opt', 'delta_n_opt', 'f_int_opt', 'tau_int_opt', 'N_opt', 'dt_opt', 'n_am_opt', 'shape_pitch_opt' ...
    , 'mu_s', 'mu_d', 'mass_obj');

%% TODO
f_fric_opt = zeros(N, 1);

for i = 1:N
    f_fric_opt(i) = -mass_obj * 9.81* (mu_d * tanh(4 * x_n_opt(2, i) / v_s) ...
       + (mu_s - mu_d) * x_n_opt(2, i) /v_s /((x_n_opt(2, i)/2/v_s)^2 + 0.75)^2 );
end
%% plot optimized trajectory
times = (1:N+1)*dt;
figure('Position',[100 400 400 400])
subplot(5,1,1)
plot(times, x_n_opt(:, :))
legend("xi", "xid")
grid on

subplot(5,1,2)
plot(times(1:N), f_fric_opt)
grid on

subplot(5,1,3)
hold on
for j = 1:n_am
    plot(times(1:N), delta_n_opt(3*j-2, :))
end
ylabel("delta x")
grid on

subplot(5,1,4)
hold on
for j = 1:n_am
    plot(times(1:N), delta_n_opt(3*j-1, :))
end
grid on
ylabel("delta y")

subplot(5,1,5)
hold on
for j = 1:n_am
    plot(times(1:N), delta_n_opt(3*j, :))
end
grid on
ylabel("delta z")

figure('Position',[500 400 800 400]);
subplot(2,3,1)
hold on
for j = 1:n_am
    plot(times(1:N), f_int_opt(3*j-2, :))
end    
grid on
title("f x")
subplot(2,3,2)
hold on
for j = 1:n_am
    plot(times(1:N), f_int_opt(3*j-1, :))
end    
grid on
title("f y")
subplot(2,3,3)
hold on
for j = 1:n_am
    plot(times(1:N), f_int_opt(3*j, :))
end    
grid on
title("f z")

subplot(2,3,4)
hold on
for j = 1:n_am
    plot(times(1:N), tau_int_opt(3*j-2, :))
end    
grid on
title("tau x")
subplot(2,3,5)
hold on
for j = 1:n_am
    plot(times(1:N), tau_int_opt(3*j-1, :))
end    
grid on
title("tau y")
subplot(2,3,6)
hold on
for j = 1:n_am
    plot(times(1:N), tau_int_opt(3*j, :))
end    
grid on
title("tau z")


%% Video
show_video = true;
save_video = false;
video_speed = 1;

if show_video
figure('Position',[500 300 600 500]);
dN = 0.1 / dt;
framesPerSecond = 1/dt/dN * video_speed;
rate = rateControl(framesPerSecond);
arrow_len = 0.2/ 9.81;


min_X = min(x_n_opt(1, :)); % Minimum for each row (x, y, z)
max_X = max(x_n_opt(1, :)); % Maximum for each row (x, y, z)

if save_video
    video_filename = '/images/test.avi';
    video = VideoWriter(video_filename);
    video.FrameRate = framesPerSecond;
    open(video);
end

for i = 1:dN:N
    clf;
    grid on; axis equal;
    set(gca, 'GridLineWidth', 1.0, 'GridAlpha', 0.3);
    blank = 1.5;
    xlim([min_X(1) - (l1+l2)/2*n_am*2.0, max_X(1) + r_o(1) + box_width * 1.5 ] )
    %ylim([min_X(3) - (l1+l2)/2*n_am*2.0, min_X(3) + (l1+l2)/2*n_am*2.0 ])
    ylim([-2.0, 2.0])
    
    
    hold on;
    R = eye(3);
    X_com_n = x_n_opt(1, i) * e_1;

    for j = 1:n_am
        R_e_j = R * R_shape{j};
        X_quad = X_com_n + R * r_cj{j};

        dx = -delta_n_opt(3*j-2, i) * arrow_len;
        dz = (-delta_n_opt(3*j, i) + mass_am(j)*9.81) * arrow_len;       
        plot(X_quad(1), X_quad(3), 'o', 'Color', 'b');
        quiver(X_quad(1), X_quad(3), dx, dz, 0, 'r', 'LineWidth', 1.5);
        
        X1 = X_quad + l1 * R_e_j * e_1;
        X2 = X_quad - l2 * R_e_j * e_1;
        plot([X1(1), X2(1)], [X1(3), X2(3)], 'b-', 'LineWidth', 1.0);
    end
    
    % com
    plot(X_com_n(1), X_com_n(3), 'o', 'Color', 'b', 'MarkerSize', 10);

    X_object = X_com_n + R * r_o;
    plot(X_object(1), X_object(3), 'o', 'Color', 'blue', 'LineWidth', 1.5);
    cornerX = box_width/2 * [-1, 1, 1, -1, -1] + X_object(1) * ones(1,5);
    cornerZ = box_height/2 * [-1, -1, 1, 1, -1] + X_object(3) * ones(1,5);
    plot(cornerX, cornerZ, 'b', 'LineWidth', 1.0);
    
    xlabel('X position'); ylabel('Z position');
    title_string = sprintf("time : %.2f sec", i*dt);
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