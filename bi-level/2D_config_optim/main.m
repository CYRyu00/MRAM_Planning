clear; close all;
addpath("../../params", "../../../casadi-3.6.7-windows64-matlab2018b")
import casadi.*
params = define_params_ver2();
mu = params{3}; r = params{4}; d = params{5};
thrust_limit = params{6}; mb = params{17}; mt = params{20}; ma = params{23};
m0 = mt + ma + mb;

B_tau = [r r -r -r;  -r r r -r; mu -mu mu -mu];
B_lambda = [1 1 1 1];

l1 = 0.35; l2 = 0.35; % ca = 0.24

e_1 = [1; 0; 0]; e_2 = [0; 1; 0]; e_3 = [0; 0; 1];

num_AMs = 2;
J_c = zeros(3, 3); % inertia w.r.t. its com
AM_com = [0; 0; 0];% 0 to com
r_0j = cell(num_AMs, 1); % 0 to j'th module
r_cj = cell(num_AMs, 1); % com to j'th module
mass_ams = m0 * ones(num_AMs, 1);
R_shape = cell(1, num_AMs);
shape_pitch = 0.0 * [10 7 3 -10 -10 10]; %[0 -10 -10 10 -10 -10 10]

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
for j = 1:num_AMs
    r_cj{j} = r_0j{j} - AM_com;% p_com to j
end

theta_min = -10 / 180 * pi * ones(num_AMs, 1);
theta_max = 10 / 180 * pi * ones(num_AMs, 1);

thrust_min = 1e-2* -thrust_limit * ones(4*num_AMs, 1);
thrust_max = thrust_limit * ones(4*num_AMs, 1);
thrust_cen = (thrust_max + thrust_min)/2; % 4 x 1
thrust_tilde = (thrust_max(1) - thrust_min(1))/2; % scalar
%% fibonacci sphere
offset = m_c* 9.81 * [0; 0.0; 1.0];
num_points = 1000; % Number of points
U = fibonacci_sphere(num_points); %  tau_y, fx, fz
beta = 1e1; % 0.01 ~ 1.0

n_feas = 10; % Number of points
target_axis = offset;
cone_angle = 1; %degree
bound_ratio = 0.7;
cone_samples = sample_3d_cone_process(target_axis, cone_angle, n_feas, bound_ratio)'; %  tau_y, fx, fz
%% Determine Q

v_major = offset / 9.81 / m_c;
r_minor = 0.1;
a = norm(v_major);
n1 = v_major / a;

N = null(n1');
V = [n1, N];
lambda1 = 1 / (a^2);
lambda2 = 1 / (r_minor^2);
lambda3 = 1 / (r_minor^2);

Lambda = diag([lambda1, lambda2, lambda3]);
Q = V * Lambda * V';
% Q = eye(3,3);
Q = Q / det(Q)^(1/length(offset));
inv_Q_half = inv(sqrtm(Q));
disp('생성된 Q^-1/2:');
disp(inv_Q_half);
% offset= [0;0;0];
%% NLP SOLVER

Theta = MX.sym('Theta', num_AMs, 1);
% center = MX.sym('center', 3, 1);
weight = MX.sym('weight', 1, 1);
f_feasible = MX.sym('f_feasible', 4*num_AMs, n_feas);
opt_variables = [Theta; weight; reshape(f_feasible, 4*num_AMs*n_feas, 1)];
opt_var0 = [(2*rand(num_AMs, 1) - 1) * 30 / 180 * pi; 1.5; repmat(thrust_cen, n_feas, 1)];

B = [];
g = [];
LBG = [];
UBG = [];
for j = 1:num_AMs
    R = R_shape{j}' * [cos(Theta(j)), 0, -sin(Theta(j)); 0, 1, 0; sin(Theta(j)), 0, cos(Theta(j))]; % R_y(-theta)
    B_j = [R, hat(r_cj{j}) * R; zeros(3, 3), R] * [B_tau; zeros(2, 4); B_lambda];
    B = [B, B_j];

    g = [g; Theta(j)];
    LBG = [LBG; theta_min(j)];
    UBG = [UBG; theta_max(j)]; 
end
B_2d = [B(2, :); B(4, :); B(6, :)];

for k = 1:n_feas
    g = [g; (cone_samples(k, :)' * norm(offset) - B_2d * f_feasible(:, k)); f_feasible(:, k)];
    LBG = [LBG; zeros(3, 1); thrust_min];
    UBG = [UBG; zeros(3, 1); thrust_max]; 
end

sum_exp = 0;
diff = B_2d * thrust_cen - offset*weight;
% diff = B_2d * thrust_cen - center;
for k = 1:num_points
    one_norm = 0;
    for i = 1:4*num_AMs
        one_norm = one_norm + abs( U(k,:)* inv_Q_half * B_2d(:, i));
        % one_norm = one_norm + sqrt((U(k,:) *inv_Q_half * B_2d(:, i))^2 + 1e-4);
    end
    support_k = U(k, :) * inv_Q_half *diff+ thrust_tilde * one_norm;
    sum_exp = sum_exp + exp(-beta * support_k);
end
obj = 1 /beta *log(sum_exp);
obj = obj + diff'*diff*5e-2;

nlp_prob = struct('x', opt_variables, 'f', obj, 'g', g);
nlp_opts = struct;
nlp_opts.ipopt.print_level = 1;
nlp_opts.ipopt.tol = 1e-5;
nlp_opts.ipopt.max_iter = 100;
nlp_opts.ipopt.mu_strategy = 'adaptive';
nlp_opts.ipopt.linear_solver = 'mumps';
nlp_opts.ipopt.jacobian_approximation = 'exact';
%nlp_opts.ipopt.hessian_approximation = 'limited-memory';
nlp_opts.print_time = 1;

solver = nlpsol('solver', 'ipopt', nlp_prob, nlp_opts);
sol = solver('x0', opt_var0, ...
             'lbx', -inf, 'ubx', inf,...
             'lbg', LBG, 'ubg', UBG);
fprintf("\nsuccess : %d", solver.stats.success)
fprintf("\nreturn_stauts : %s\n", solver.stats.return_status)
%%
solution = full(sol.x);
opt_f = full(sol.f);
theta_opt = solution(1:num_AMs);
center_opt = offset * solution(num_AMs+1);
% center_opt = solution(num_AMs+1:num_AMs+3);
fprintf("theta opt(degree) :");
disp(theta_opt' * 180 / pi)
fprintf("Objective function :");
disp(opt_f)

B_opt = [];
for j = 1:num_AMs
    R = R_shape{j}' * [cos(theta_opt(j)), 0, -sin(theta_opt(j)); 0, 1, 0; sin(theta_opt(j)), 0, cos(theta_opt(j))]; % R_y(-theta)
    B_j = [R, hat(r_cj{j}) * R; zeros(3, 3), R] * [B_tau; zeros(2, 4); B_lambda];
    B_opt = [B_opt, B_j];
end
fprintf("\ndiff :\n")
disp([B_opt(2,:); B_opt(4,:); B_opt(6,:)] * thrust_cen - center_opt)
%% plot wrench polytope
B_2d_opt = [];
for j = 1:num_AMs
    B_j_opt = 2 * [B_opt(2,4*j-3:4*j-2);
                   B_opt(4,4*j-3:4*j-2);
                   B_opt(6,4*j-3:4*j-2)];
    B_2d_opt = [B_2d_opt, B_j_opt]; %tau_y, f_x, f_z
end

n_actuators = 2*num_AMs; % 입력 개수
f_min = thrust_min(1);
f_max = thrust_max(1);

%(Vertex Enumeration)
num_vertices = 2^n_actuators;
perms = dec2bin(0:num_vertices-1) - '0';
f_vertices = perms' * (f_max - f_min) + f_min;
w_points = B_2d_opt * f_vertices; 

tau_y_arr = w_points(1, :)';
f_x_arr = w_points(2, :)';
f_z_arr = w_points(3, :)';

[K, volume] = convhull(tau_y_arr, f_x_arr, f_z_arr);

figure('Color', 'w');
trisurf(K, tau_y_arr, f_x_arr, f_z_arr, 'FaceColor', 'cyan', 'FaceAlpha', 0.3, 'EdgeColor', 'b');
hold on;

grid on; axis equal;
xlabel('Torque (\tau_y)');
ylabel('Force X (f_x)');
zlabel('Force z (f_z)');
title(['Wrench Polytope (Vol: ', num2str(volume, '%.2f'), ')']);

plot3(0,0,0, 'k+', 'MarkerSize', 15, 'LineWidth', 2);
plot3(center_opt(1), center_opt(2), center_opt(3), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
plot3(offset(1), offset(2), offset(3), 'r+', 'MarkerSize', 15, 'LineWidth', 2);



[x_s, y_s, z_s] = sphere(30);
sphere_coords = [x_s(:)'; y_s(:)'; z_s(:)'];

ellipsoid_coords = inv_Q_half * sphere_coords * abs(opt_f) + center_opt;

xe = reshape(ellipsoid_coords(1,:), size(x_s));
ye = reshape(ellipsoid_coords(2,:), size(y_s));
ze = reshape(ellipsoid_coords(3,:), size(z_s));

surf(xe, ye, ze, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

axis equal; grid on;