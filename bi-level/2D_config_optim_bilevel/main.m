clear; close all
addpath("../../params")

warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:illConditionedMatrix');
%% parameters
num_points = 100; % Number of points
num_AMs = 5;

params = define_params_ver2();
mu = params{3}; r = params{4}; d = params{5};
thrust_limit = params{6}; mb = params{17}; mt = params{20}; ma = params{23};
m0 = mt + ma + mb;

B_tau = [r r -r -r;  -r r r -r; mu -mu mu -mu];
B_lambda = [1 1 1 1];

l1 = 0.35; l2 = 0.35;

e_1 = [1; 0; 0]; e_2 = [0; 1; 0]; e_3 = [0; 0; 1];

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

thrust_min = 1e-2*-thrust_limit; %thrusts limit
thrust_max = thrust_limit;
%% optimization variables
% theta0 = (2*rand(num_AMs, 1) - 1) * 5 / 180 * pi;
theta0 = [1; 2; 3; 4; 5]/ 180 * pi;

lb = -70 / 180 * pi * ones(num_AMs, 1);
ub = 70 / 180 * pi * ones(num_AMs, 1);

offset = [0; 0.6; 1]; % ty, fx, fz
center_axis = offset;
angle_deg = 5;
num_samples = num_points;
bound_ratio = 0.7;
U = sample_3d_cone_process(center_axis, angle_deg, num_samples, bound_ratio);
%% gradient check 
% checkGradients(@(theta)test_const(theta, num_AMs, f_min, f_max), theta0, Display="on", IsConstraint=true);
%clc
% [valid, err]=checkGradients(@(theta)high_cost(theta, U, num_AMs, num_points, thrust_min, thrust_max, R_shape, r_cj, B_tau, B_lambda), theta0, Display="on", IsConstraint=false);
%% 
options = optimoptions('fmincon', ...
    'SpecifyObjectiveGradient', true, ...
    'UseParallel', true, ...
    'Display', 'iter-detailed', ...
    'Algorithm', 'interior-point', ... % active-set, interior-point
    'HessianApproximation', 'lbfgs', ...
    'OptimalityTolerance', 1e-6, ...
    'ConstraintTolerance', 1e-6, ...
    'StepTolerance', 1e-6);
% tic
% [thetasol, fval, exflag, output, lambda, grad, hessian] = ...
%   fmincon( @(theta)test_cost(theta, num_AMs, num_points), theta0, [], [], [], [], lb, ub, [], options);
% toc

tic
[theta_opt, fval, exflag, output, lambda, grad, hessian] = ...
  fmincon( @(theta)high_cost(theta, U, num_AMs, num_points, thrust_min, thrust_max, R_shape, r_cj, B_tau, B_lambda), theta0, [], [], [], [], lb, ub, [], options);
toc

fprintf("\ntheta_opt : ")
disp(theta_opt' * 180/pi)
%% plot wrench polytope
B_opt = [];
for j = 1:num_AMs
    R = R_shape{j}' * [cos(theta_opt(j)), 0, -sin(theta_opt(j)); 0, 1, 0; sin(theta_opt(j)), 0, cos(theta_opt(j))]; % R_y(-theta)
    B_j = [R, hat(r_cj{j}) * R; zeros(3, 3), R] * [B_tau; zeros(2, 4); B_lambda];
    B_opt = [B_opt, B_j];
end

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
plot3(offset(1), offset(2), offset(3), 'r+', 'MarkerSize', 15, 'LineWidth', 2);


[sx, sy, sz] = sphere(50);
surf(sx, sy, sz, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

c_norm = center_axis / norm(center_axis)* abs(fval);
quiver3(0, 0, 0, c_norm(1)*1.2, c_norm(2)*1.2, c_norm(3)*1.2, ...
    'g', 'LineWidth', 1, 'MaxHeadSize', 0.2, 'DisplayName', 'Target Axis');

scatter3(U(1,:) * abs(fval), U(2,:)* abs(fval), U(3,:)* abs(fval), 15, 'r', 'filled', ...
    'MarkerFaceAlpha', 0.3, 'DisplayName', 'Samples');

axis equal; grid on;