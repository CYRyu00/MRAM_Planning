function [theta_opt, does_success, opt_f, Q]  = get_opt_config(mean, cov,  r_min, r_max, theta_init, theta_min, theta_max, num_points, beta, theta, theta_tilde, thrust_min, thrust_max, m0, mass_ams, R_shape, l1, l2, B_tau, B_lambda, num_AMs, r_cj)
addpath("../../../casadi-3.6.7-windows64-matlab2018b")
import casadi.*

e_1 = [1; 0; 0]; e_2 = [0; 1; 0]; e_3 = [0; 0; 1];

theta_min = max(theta_min, theta - theta_tilde * ones(num_AMs, 1));
theta_max = min(theta_max, theta + theta_tilde * ones(num_AMs, 1));
% disp(theta /pi * 180)
% disp(theta_min /pi * 180)
% disp(theta_max /pi * 180)
thrust_cen = (thrust_max + thrust_min)/2; % 4 x 1
thrust_tilde = (thrust_max(1) - thrust_min(1))/2; % scalar

%% fibonacci sphere
offset = mean;
U = fibonacci_sphere(num_points); %  tau_y, fx, fz
k_diff = 1e-1 * num_AMs;

%% Determine Q
v_major = offset;
a = norm(v_major);
n1 = v_major / a;

N = null(n1');
V = [n1, N];
cov2 = N(:, 1)' * cov * N(:, 1);
cov3 = N(:, 2)' * cov * N(:, 2);
k_cov = 1e-2;
lambda1 = 1 / (1^2);
lambda2 = 1 / min(r_max^2, max(r_min^2, k_cov * cov2));
lambda3 = 1 / min(r_max^2, max(r_min^2, k_cov * cov3));

Lambda = diag([lambda1, lambda2, lambda3]);
Q = V * Lambda * V';
Q = Q / det(Q)^(1/length(offset));
% Q = eye(3,3);

inv_Q_half = inv(sqrtm(Q));
% disp('생성된 Q^-1/2:');
% disp(inv_Q_half);
% offset= [0;0;0];
%% NLP SOLVER
Theta = MX.sym('Theta', num_AMs, 1);
f_feas = MX.sym('f_feas', 4*num_AMs, 1);
weight = MX.sym('weight', 1, 1);
opt_variables = [Theta; f_feas; weight];
opt_var0 = [theta_init; thrust_cen; 1];

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

g = [g; (offset - B_2d * f_feas); f_feas];
LBG = [LBG; zeros(3, 1); thrust_min];
UBG = [UBG; zeros(3, 1); thrust_max]; 

sum_exp = 0;
for k = 1:num_points
    one_norm = 0;
    for i = 1:4*num_AMs
        one_norm = one_norm + abs( U(k,:)* inv_Q_half * B_2d(:, i));
        % one_norm = one_norm + sqrt((U(k,:) *inv_Q_half * B_2d(:, i))^2 + 1e-4);
    end
    support_k = thrust_tilde * one_norm;
    sum_exp = sum_exp + exp(-beta * support_k);
end
obj = 1 /beta *log(sum_exp);

diff = B_2d*thrust_cen - weight*offset;
obj = obj + diff' * diff * k_diff; 

nlp_prob = struct('x', opt_variables, 'f', obj, 'g', g);
nlp_opts = struct;
nlp_opts.ipopt.print_level = 1;
nlp_opts.ipopt.tol = 1e-5;
nlp_opts.ipopt.max_iter = 100;
nlp_opts.ipopt.mu_strategy = 'adaptive';
nlp_opts.ipopt.linear_solver = 'mumps';
nlp_opts.ipopt.jacobian_approximation = 'exact';
%nlp_opts.ipopt.hessian_approximation = 'limited-memory';
nlp_opts.print_time = 0;

solver = nlpsol('solver', 'ipopt', nlp_prob, nlp_opts);
sol = solver('x0', opt_var0, ...
             'lbx', -inf, 'ubx', inf,...
             'lbg', LBG, 'ubg', UBG);
% fprintf("\nsuccess : %d", solver.stats.success)
does_success = solver.stats.success;
% fprintf("\nreturn_stauts : %s\n", solver.stats.return_status)

%%
solution = full(sol.x);
opt_f = full(sol.f);
theta_opt = solution(1:num_AMs);

% B_opt = [];
% for j = 1:num_AMs
%     R = R_shape{j}' * [cos(theta_opt(j)), 0, -sin(theta_opt(j)); 0, 1, 0; sin(theta_opt(j)), 0, cos(theta_opt(j))]; % R_y(-theta)
%     B_j = [R, hat(r_cj{j}) * R; zeros(3, 3), R] * [B_tau; zeros(2, 4); B_lambda];
%     B_opt = [B_opt, B_j];
% end
% center_opt = [B_opt(2,:); B_opt(4,:); B_opt(6,:)] * thrust_cen;
% diff_opt = center_opt - solution(end) * offset;
% radius_opt = - opt_f + diff_opt'*diff_opt * k_diff;

% fprintf("\ntheta opt(degree) :");
% disp(theta_opt' * 180 / pi)
% fprintf("radius: %.4f\n", radius_opt);

end