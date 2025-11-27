clear; close all;
addpath("../params",  "../../../casadi-3.6.7-windows64-matlab2018b")
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
l1 = 0.4; l2 = 0.3;
kappa = 30;
%%
M_am = zeros(n_am * 6, n_am * 6);
for j =1:n_am
    M_am(3*j-2: 3*j, 3*j-2: 3*j) = mass_am(j) * eye(3);
    M_am(3*(j+n_am)-2: 3*(j+n_am), 3*(j+n_am)-2: 3*(j+n_am)) = It + Ib;
end
M_o = diag(mass_obj * ones(3,1));
M = [M_o, zeros(3, n_am*6); zeros(n_am*6, 3), M_am];

p = [ones(3, 1) * 0; ones(n_am*3, 1)];
phi = ones(n_am, 1)*0;
qd = ones(3 + n_am*6, 1) * 0;
U = [ones((n_am)*3, 1) * 10; ones(n_am*3, 1) * 0];

[qdd, f_in, tau_in] = FD_cart(p, phi, qd, U, M, k_obj, b_obj, l1, l2, n_o, n_am);
%% Dynamics
p = MX.sym('p', n_o + n_am*3, 1); % [p_o; p_am];
phi = MX.sym('phi', n_am, 1);
qd = MX.sym('qd', n_o + n_am*6, 1); % [pd_o; pd_am; w_am]
u = MX.sym('u', n_am*3, 1); % [force_am]
tau = Mx.zeros(n_am*3,1); % tau_am

U = [u; tau];
[qdd, f_in, tau_in] = FD_cart(p, phi, qd, U, M, k_obj, b_obj, l1, l2, n_o, n_am);

qd_next = qd + qdd * dt;
p_next = p + qd(1:n_o + n_am*3, 1) *dt;
for j = 1:n_am  
    phi_next(j) = phi(j) + [0, 1, 0] * qd(3*j-2:3*j) *dt;
end

fd = Function('fd', {p, phi, qd, u}, {p_next, phi_next, qd_next});]


Delta_n = MX.sym('Delta_n', n_am*3, N);
Q_n = MX.sym('Q_n', n_o*2 + n_am*10, N + 1); %[p;phi;qd]
opt_variables = [reshape(X, nx * (N + 1), 1); reshape(U, n_am*3 * N, 1)];


g = [sum(Delta_p, 2); moment];

opt_var = [Delta_p(:)]; 

nlp_prob = struct('x', opt_var, 'f', obj, 'g', g);
nlp_opts = struct;
nlp_opts.ipopt.print_level = 0;
nlp_opts.ipopt.tol = 1e-3;
%nlp_opts.ipopt.max_iter = max_iter;
nlp_opts.ipopt.mu_strategy = 'adaptive';

nlp_opts.ipopt.jacobian_approximation = 'exact';
% nlp_opts.ipopt.hessian_approximation = 'limited-memory';
nlp_opts.print_time = 0;

solver = nlpsol('solver', 'ipopt', nlp_prob, nlp_opts);
opt_var0 = [Delta_p_opt(:)];
LBG = delta_hat_tot + [0;0;0; e_2' * (k_R_cen * e_R_cen + k_w_cen * e_w_cen)];
UBG = LBG;
LBG = [LBG; -inf * ones(num_AMs*2, 1)];
UBG = [UBG; 0 * ones(num_AMs*2, 1)];
sol = solver('x0', opt_var0, 'lbg', LBG, 'ubg', UBG);