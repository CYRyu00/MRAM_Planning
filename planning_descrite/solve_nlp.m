function [x_opt, u_opt, optimal_value, exit_flag, processing_time] ....
    =  solve_nlp(params, shape, num_AMs, mass, inertia, r_i_ci, q_o_ref, tau_min, tau_max, u_min, u_max, x_0, x_f, X_init_guess, dt, N)
import casadi.*

m0 = params{1}; I0 = params{2}; mu = params{3}; r= params{4}; d= params{5};
thrust_limit= params{6}; kt = params{7}; c_1 = params{8}; c_2 = params{9}; mass_door = params{10};
handle_factor = params{11}; n = params{14}; dh = params{15}; gravity = params{16};

nx = n * 2;
nu = 2 + num_AMs * 4;
%% Dynamics
x = MX.sym('x', nx, 1); % [q; qd];
u = MX.sym('u', nu, 1);

tau = [(-c_1 * x(5)); (-c_2 * x(6) -kt * x(2) + mass{2} * handle_factor); u(1); u(2)] ;
F_ext = map_u2wrench( u(3:end), shape, mu, r, d);

qdd = FD_ver2(n, dh, mass, inertia, r_i_ci, gravity, x(1:4), x(5:8), tau, F_ext);

x_dot = [x(5:8); qdd];
x_next = x + dt * x_dot ;
f = Function('f', {x, u}, {x_next});
U = MX.sym('U', nu, N);
X = MX.sym('X', nx, N+1);

opt_variables = [reshape(X, nx * (N + 1), 1); reshape(U, nu * N, 1)];
%%
obj = 0;
g = [];
[Q1, Q2, Q3, R] = define_nlp_params(nu);

for k = 1:N
    if k <= 1 / dt
        obj = obj + (X(:, k) - x_0)' * Q3 * (X(:, k) - x_0);
    elseif N-k <= 1 / dt
        obj = obj + (X(:, k) - x_f)' * Q3 * (X(:, k) - x_f);
    else
        obj = obj + X(:, k)' * Q1 * X(:, k) + (X(1:4, k) - q_o_ref(:, k))' * Q2 * (X(1:4, k) - q_o_ref(:, k));
        obj = obj + U(:, k)' * R * U(:, k);
    end
    % dynamics equality constraint
    g = [g; X(:, k + 1) - f(X(:, k), U(:, k))];
end
g = [g; X(:, 1) - x_0];
g = [g; X(:, N + 1) - x_f];

%%
nlp_prob = struct('x', opt_variables, 'f', obj, 'g', g);
nlp_opts = struct;
nlp_opts.ipopt.print_level = 1;
nlp_opts.ipopt.tol = 1e-3;
nlp_opts.ipopt.max_iter = 2000;
nlp_opts.ipopt.mu_strategy = 'adaptive';
nlp_opts.ipopt.linear_solver = 'mumps';
%nlp_opts.ipopt.jacobian_approximation = 'exact';
nlp_opts.ipopt.hessian_approximation = 'limited-memory';
nlp_opts.print_time = 1;

solver = nlpsol('solver', 'ipopt', nlp_prob, nlp_opts);
% solve
LBx = [repmat([ - 0.1 * pi; - pi/3; - pi/2 ; - pi/2; - inf; - inf; - inf; - inf], N + 1, 1); ...
       repmat([tau_min; tau_min; u_min*ones(nu-2, 1)], N, 1)];
UBx = [repmat([ pi/2; pi/3; pi/2 ; pi/2; inf; inf; inf; inf], N + 1, 1); ...
       repmat([tau_max; tau_max; u_max * ones(nu - 2,1)], N, 1)];

sol = solver('x0',  X_init_guess, ...
    'lbx', LBx, 'ubx', UBx,...
    'lbg', 0, 'ubg', 0);

solution = full(sol.x);
x_opt = reshape(solution(1:nx * (N + 1)), nx, N + 1)';
u_opt = reshape(solution(nx*(N + 1) + 1:end), nu, N)';
optimal_value = full(sol.f);
exit_flag = solver.stats.success;
processing_time = solver.stats.t_wall_total;
disp(solver.stats.return_status)
end