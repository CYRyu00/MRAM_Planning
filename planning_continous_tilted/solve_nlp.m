function [rho_opt, x_opt, u_opt, optimal_value,exit_flag,processing_time, rho_history, exit_flag_history, time_history] ....
    =  solve_nlp(params, theta, num_AMs, K, L, core, q_o_ref, tau_min, tau_max, u_min, u_max, x_0, x_f, X_init_guess, dt, N, t1, max_iter, eps, gamma, Q1, Q2, Q3, R)
%% NLP SOLVER
import casadi.*

m0 = params{1}; I0 = params{2}; mu = params{3}; r = params{4}; d = params{5};
thrust_limit= params{6}; kt = params{7}; c_1 = params{8}; c_2 = params{9}; mass_door = params{10};
handle_factor = params{11}; inertia = params{12}; r_i_ci = params{13}; n = params{14}; dh = params{15}; gravity = params{16};

nx = n * 2;
nu = K * L * 8;

x = MX.sym('x', nx, 1); % [q; qd];
u = MX.sym('u', nu, 1);
rho = MX.sym('rho', K, L); %shape

%% Dynamics
[AM_com, AM_mass, AM_inertia] = get_inertia_duo(rho, K, L, core, m0, I0, d);
mass = {mass_door(1), mass_door(2), mass_door(3), AM_mass};
inertia{4} = AM_inertia;
r_i_ci{4} = [AM_com(1); r_i_ci{4}(2); AM_com(2)];

tau = [(- c_1 * x(5)); (- c_2 * x(6) - kt * x(2) + mass{2} * handle_factor); 0; 0] ;
F_ext = map_u2wrench(u, rho, K, L, core, mu, r, d, theta);

qdd = FD(n, dh, mass, inertia, r_i_ci, gravity, x(1:4), x(5:8), tau, F_ext);

x_dot = [x(5:8); qdd];
x_next = x + dt * x_dot ;

f = Function('f', {x, u, rho}, {x_next});

U = MX.sym('U', nu, N);
X = MX.sym('X', nx, N + 1);
opt_variables = [reshape(rho, K * L,1); reshape(X, nx * (N + 1), 1); reshape(U, nu * N, 1)];
%% Define NLP and solve iteratively
success = 0;
iter = 1;

rho_history = cell(1,100);
exit_flag_history = cell(1,100);
time_history = cell(1,100);
obj_history = cell(1,100);

tic
while( success==0 || eps> (0.005 * gamma)) && iter < 15
    fprintf("num AMs = %d, iter = %d, eps = %.5f\n", num_AMs, iter, eps);
    obj = 0;
    g = [];

    for k = 1:N
        if k <= t1 / dt
            obj = obj + (X(:, k) - x_0)' * Q3 * (X(:, k) - x_0);
        elseif N - k <= t1 / dt
            obj = obj + (X(:, k) - x_f)' * Q3 * (X(:, k) - x_f);
        else
            obj = obj + X(:, k)' * Q1 * X(:, k) + (X(1:4, k) - q_o_ref(:, k))' * Q2 * (X(1:4, k) - q_o_ref(:, k));
            obj = obj + U(:, k)' * R * U(:, k);
        end
        % dynamics equality constraint
        g = [g; X(:, k + 1) - f(X(:, k), U(:, k), rho) ; - X(:, k + 1) + f(X(:, k), U(:, k), rho)];
    end

    for i_=1:K
        for j_=1:L
            if(i_ == core(1) && j_ == core(2))
                g = [g; rho(i_, j_) - 1; - rho(i_, j_) + 1];
            else
                neighbor_max = MX(0);
                if i_ > 1 && i_>core(1)
                    neighbor_max = fmax(neighbor_max, rho(i_ - 1, j_));
                end
                if i_ < K && i_<core(1)
                    neighbor_max = fmax(neighbor_max, rho(i_ + 1, j_));
                end
                if j_ > 1 && j_>core(2)
                    neighbor_max = fmax(neighbor_max, rho(i_, j_ - 1));
                end
                if j_ < L && j_<core(2)
                    neighbor_max = fmax(neighbor_max, rho(i_, j_ + 1));
                end
                g = [g; - rho(i_, j_); rho(i_, j_) - neighbor_max ];
            end
            g = [g; rho(i_,j_) * (rho(i_,j_) - 1) - eps; - rho(i_,j_) * (rho(i_,j_) - 1) - eps];
        end
    end
    g = [g; sum1(sum2(rho)) - num_AMs; - sum1(sum2(rho)) + num_AMs ];

    g = [g; X(:, 1) - x_0; - X(:, 1) + x_0];
    g = [g; X(:, N + 1) - x_f ; -X(:, N + 1) + x_f];

    % solver options
    nlp_prob = struct('x', opt_variables, 'f', obj, 'g', g);
    nlp_opts = struct;
    nlp_opts.ipopt.print_level = 1;
    nlp_opts.ipopt.tol = 1e-3;
    nlp_opts.ipopt.max_iter = max_iter;
    nlp_opts.ipopt.mu_strategy = 'monotone';
    nlp_opts.ipopt.linear_solver = 'mumps';
    %nlp_opts.ipopt.jacobian_approximation = 'exact';
    nlp_opts.ipopt.hessian_approximation = 'limited-memory';
    nlp_opts.print_time = 1;

    solver = nlpsol('solver', 'ipopt', nlp_prob, nlp_opts);

    % solve
    LBx = [ ones(K * L, 1) * - inf; ...
        repmat([ - 0.1 * pi; - pi/3; - pi/2 ; - pi/2; - inf; - inf; - inf; - inf], N + 1, 1); ...
        repmat(u_min * ones(nu, 1), N, 1)];
    UBx = [ ones(K * L, 1) * inf; ...
        repmat([ pi/2; pi/3; pi/2 ; pi/2; inf; inf; inf; inf], N + 1, 1); ...
        repmat(u_max * ones(nu,1), N, 1)];

    if iter==1
        sol = solver('x0',  X_init_guess, ...
            'lbx', LBx, 'ubx', UBx,...
            'lbg', -inf, 'ubg', 0);
    else
        sol = solver('x0',  sol.x, ...
            'lbx', LBx, 'ubx', UBx,...
            'lbg', -inf, 'ubg', 0);
    end

    full_sol = full(sol.x);
    current_rho = reshape(full_sol(1:K * L), K, L);
    rho_history{iter} = current_rho;
    exit_flag_history{iter} = solver.stats.return_status;
    time_history{iter} = solver.stats.t_wall_total;
    obj_history{iter} = full(sol.f);
    fprintf("current object function value: %f\n", full(sol.f));
    disp(solver.stats.return_status)
    disp(current_rho)
    
    success = solver.stats.success;
    eps = eps * gamma;
    iter= iter + 1;
end
toc

solution = full(sol.x);
rho_opt =reshape(solution(1:K * L), K, L);
x_opt = reshape(solution(K * L + 1 : K * L + nx * (N + 1)), nx, N + 1)';
u_opt = reshape(solution((K * L + nx * (N + 1) + 1 ):end), nu, N)';
optimal_value = full(sol.f);
exit_flag = solver.stats.success;
processing_time = sum(toc);
disp(solver.stats.return_status)
end
