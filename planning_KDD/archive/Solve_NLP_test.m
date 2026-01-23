clear; close all;
addpath("../params",  "../../casadi-3.6.7-windows64-matlab2018b")
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

mass_obj = 10;
k_obj = 20;
b_obj = 30; % crt:28.28
n_o = 3;

dt = 0.1; N = 10/dt;
n_am = 4;
mass_am = 2.0 * ones(1,n_am); m_com = sum(mass_am);
l1 = 0.35; l2 = 0.3;

p_1_i = [0;0;0; 0;0;0];
p_1_f = [-0.3;0;0; 0;0;0];
tan_max = tan(45 / 180 *pi);
f_int_max = 100;
tau_int_max = 100;
k_smooth = 1e-1;

max_iter = 10000;
%%
M_am = zeros(n_am * 6, n_am * 6);
for j =1:n_am
    M_am(3*j-2: 3*j, 3*j-2: 3*j) = mass_am(j) * eye(3);
    M_am(3*(j+n_am)-2: 3*(j+n_am), 3*(j+n_am)-2: 3*(j+n_am)) = It + Ib;
end
M_o = diag(mass_obj * ones(3,1));
M = [M_o, zeros(3, n_am*6); zeros(n_am*6, 3), M_am];

shape_pitch =[0 0 0 0 0 0] / 180 *pi; % = [0 -10 -20 -20 -10 0 10] / 180 *pi;
r_g = [-0.3; 0; 0.2];

%% Dynamics
p_1 = MX.sym('p_1', 3); % [p_1];
pd_1 = MX.sym('pd_1', 3);
phi_1 = MX.sym('phi_1', 1);
w_1 = MX.sym('phi_1', 3);
u = MX.sym('u', n_am*3, 1); % gravity compensated force lamada R e3 = u + mge3, u = -delta
tau = MX.zeros(n_am*3,1); % tau_am

U = [u; tau];
[qdd, f_int, tau_int] = FD_cart_casadi(p_1, phi_1, pd_1, w_1, shape_pitch, r_g, U, M, k_obj, b_obj, l1, l2, n_o, n_am);

pd_1_next = pd_1 + qdd(4:6) *dt;
p_1_next = p_1 + pd_1_next *dt;
w_1_next = w_1 + qdd(4 + 3*n_am : 6 + 3*n_am);
phi_1_next = phi_1 + [0, 1, 0] * w_1_next * dt;

fd = Function('fd', {p_1, pd_1, phi_1, w_1, u}, {p_1_next, pd_1_next, phi_1_next, w_1_next});
f_int_func = Function('f_int_func', {p_1, pd_1, phi_1, w_1, u}, {f_int});
tau_int_func = Function('tau_int_func', {p_1, pd_1, phi_1, w_1, u}, {tau_int});

delta_n = MX.sym('delta_n', n_am * 3, N); % norminal delta = -u
p_n = MX.sym('p_1_n', 6, N + 1); % norminal trajectory [p_1; pd_1]
phi_n = MX.sym('phi_n', 1, N + 1);
w_n = MX.sym('w_n', 3, N + 1);
t = MX.sym('t', 2, N ); % inf norm of 2-norm of thrusts, smoothing
opt_variables = [reshape(t, 2 * N, 1); reshape(p_n, 6 * (N + 1), 1); reshape(phi_n, N + 1, 1); ...
                reshape(w_n, 3 * (N + 1), 1); reshape(delta_n, n_am * 3 * N, 1)];
opt_var0 = [zeros(2 * N, 1); zeros(6 * (N + 1), 1); zeros(N + 1, 1); ...
            zeros(3 * (N + 1), 1); zeros(n_am * 3 * N, 1)];
%% Define NLP

% start position, end position
g = [p_n(:, 1); p_n(:, end)];
LBG = [p_1_i; p_1_f]; 
UBG = [p_1_i; p_1_f];

g = [g; delta_n(:, 1)];
LBG = [LBG; zeros(3*n_am, 1)];
UBG = [UBG; zeros(3*n_am, 1)];

obj = 0;
for i = 1: N
    % %tiling angle ineq constraints
    % for j=1:n_am
    %     Rj = Ry(phi_n(i) + shape_pitch(j));
    %     g = [g; (Rj*e_1)' * (mass_am(j) * 9.81 * e_3 - delta_n(3*j-2:3*j, i)) ];
    %     LBG = [LBG; -tan_max * (Rj*e_3)' * (mass_am(j) * 9.81 * e_3 - delta_n(3*j-2:3*j, i))];
    %     UBG = [UBG; tan_max * (Rj*e_3)' * (mass_am(j) * 9.81 * e_3 - delta_n(3*j-2:3*j, i))];
    % end
    % %internal wrench constraints
    % g = [g; f_int_func(p_n(1:3, i), p_n(4:6, i), phi_n(i), w_n(:,i), -delta_n(:, i)); ...
    %         tau_int_func(p_n(1:3, i), p_n(4:6, i), phi_n(i), w_n(:,i), -delta_n(:, i))];
    % LBG = [LBG; -f_int_max * ones(3*n_am, 1); -tau_int_max * ones(3*n_am, 1)]; 
    % UBG = [UBG; f_int_max * ones(3*n_am, 1); tau_int_max * ones(3*n_am, 1)];
    
    % Dynamics inequality
    [p_1_next, pd_1_next, phi_1_next, w_1_next] = fd(p_n(1:3, i), p_n(4:6, i), phi_n(i), w_n(:,i), -delta_n(:, i));
    g = [g; [p_n(:, i+1); phi_n(i+1); w_n(:,i+1)] - [p_1_next; pd_1_next; phi_1_next; w_1_next] ];
    LBG = [LBG; zeros(10, 1)];
    UBG = [UBG; zeros(10, 1)];

    % object funtion
    % for j = 1:n_am
    %     feedforward = mass_am(j) * 9.81 * e_3 - delta_n(3*j-2:3*j, i);
    %     obj = obj + feedforward' * feedforward;
    % end
    
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
    %     g = [g; k_smooth * (delta_n(:, i)- delta_n(:, i-1))];
    %     LBG = [LBG; -t(2, i) * ones(3*n_am,1)];
    %     UBG = [UBG; t(2, i) * ones(3*n_am,1)];
    % end
end

% solver options
nlp_prob = struct('x', opt_variables, 'f', obj, 'g', g);
nlp_opts = struct;
nlp_opts.ipopt.print_level = 5;
nlp_opts.ipopt.tol = 1e-3;
nlp_opts.ipopt.max_iter = max_iter;
nlp_opts.ipopt.mu_strategy = 'monotone';
nlp_opts.ipopt.linear_solver = 'mumps';
nlp_opts.ipopt.jacobian_approximation = 'exact';
% nlp_opts.ipopt.hessian_approximation = 'limited-memory';
nlp_opts.print_time = 5;

solver = nlpsol('solver', 'ipopt', nlp_prob, nlp_opts);
sol = solver('x0', opt_var0, ...
             'lbx', -inf, 'ubx', inf,...
             'lbg', LBG, 'ubg', UBG);
full(sol.f)

solution = full(sol.x);
t_opt = reshape(solution(1:2*N), 2, N); % inf norm of 2-norm of thrusts, smoothing
p_n_opt = reshape(solution(2*N+1: 2*N + 6*(N + 1)), 6, N+1); % norminal trajectory [p_1; pd_1]
phi_n_opt = reshape(solution(2*N + 6*(N + 1) + 1: 2*N + 7*(N + 1)), 1, N+1);
w_n_opt = reshape(solution(2*N + 7*(N + 1) + 1: 2*N + 10*(N + 1)), 3, N+1);
delta_n_opt = reshape(solution(2*N + 10*(N + 1) + 1: end), 3*n_am, N); % norminal delta = -u

opt_variables = [reshape(t, 2 * N, 1); reshape(p_n, 6 * (N + 1), 1); reshape(phi_n, N + 1, 1); ...
                reshape(w_n, 3 * (N + 1), 1); reshape(delta_n, n_am * 3 * N, 1)];

disp(solver.stats.success)
function out = Ry(q)
    out = [cos(q), 0, sin(q);
           0, 1, 0;
           -sin(q), 0, cos(q)];
end