clear; close all;
addpath("../functions", "../../params")

%% Dynamics Parameters
params = define_params_ver2();
mu = params{3}; r = params{4};
thrust_limit = params{6};
B = [r r -r -r;  -r r r -r; mu -mu mu -mu; 1 1 1 1];

n_am = 4;
mass = 2.0 * ones(1,n_am); m_com = sum(mass);
l1 = 0.4; l2 = 0.3;

e_1 = [1;0;0]; e_2 = [0;1;0]; e_3 = [0;0;1];

a = 0.2; b = 0.2;c = 0.1;
J_quad = diag([b^2+c^2, a^2+c^2, a^2+b^2]) /5 * mass(1) * 0.7;

a = (l1+l2)/2; b = 0.05; c = 0.05;
J_tool = diag([b^2+c^2, a^2+c^2, a^2+b^2]) /5 * mass(1) * 0.3;

%% 
M = zeros(n_am * 6, n_am * 6);
for j =1:n_am
    M(3*j-2: 3*j, 3*j-2: 3*j) = mass(j) * eye(3);
    M(3*(j+n_am)-2: 3*(j+n_am), 3*(j+n_am)-2: 3*(j+n_am)) = J_tool + J_quad;
end

A = zeros((n_am-1)*6,n_am*6);
R1 = eye(3);
R2 = eye(3); % TODO

for j = 1:n_am - 1
    A(3*j-2: 3*j, 3*j-2: 3*j+3) = [eye(3), -eye(3)];
    A(3*j-2: 3*j, 3*(j+n_am)-2: 3*(j+n_am)+3) = [l2 * R1 * S(e_1), l1 * R2 * S(e_1)];
    A(3*(j+n_am-1)-2: 3*(j+n_am-1), 3*(j+n_am)-2: 3*(j+n_am)+3) = [eye(3), -eye(3)];
end

% A_dagger = - A' * inv(A* inv(M) * A') * A * inv(M)
A_dagger = - A' * ((A* (M\ A')) \ A)* inv(M);
A_dagger_2 = ((A* (M\ A')) \ A)* inv(M);
A_dagger_3 = (A* A') \ A;

% (A_dagger_2 - A_dagger_3) * [ones(n_am*3, 1) * 1; ones(n_am*3,1)* 1]
%% Internal forces 

A_f_t = zeros(n_am * 3, (n_am - 1) * 3);
A_f_r = zeros(n_am * 3, (n_am - 1) * 3);
A_tau_r = zeros(n_am * 3, (n_am - 1) * 3);
for j = 1 : n_am % f2, f3, f4 ...
    Rj = eye(3);
    if j == 1
        A_f_t(3*j-2:3*j, 3*j-2:3*j) = -eye(3);
        A_f_r(3*j-2:3*j, 3*j-2:3*j) = -l2*S(Rj*e_1);
        A_tau_r(3*j-2:3*j, 3*j-2:3*j) = -eye(3);
    elseif j == n_am
        A_f_t(3*j-2:3*j, 3*j-5:3*j-3) = eye(3);
        A_f_r(3*j-2:3*j, 3*j-5:3*j-3) = -l1*S(Rj*e_1);
        A_tau_r(3*j-2:3*j, 3*j-5:3*j-3) = eye(3);
    else
        A_f_t(3*j-2:3*j, 3*j-5:3*j) = [eye(3), -eye(3)];
        A_f_r(3*j-2:3*j, 3*j-5:3*j) = [-l1*S(Rj*e_1), -l2*S(Rj*e_1)];
        A_tau_r(3*j-2:3*j, 3*j-5:3*j) = [eye(3), -eye(3)];
    end
end
%
A_u2f = A_f_t(4:end, :)\A_dagger(4 : n_am*3, :);
A_u2tau = A_tau_r(4:end, :)\A_f_r(4:end, :)*A_u2f + A_tau_r(4:end, :)\ A_dagger((n_am+1)*3+1:end, :);
%%
input = [ones(n_am*3, 1) * 1; ones(n_am*3,1)* 1];
force = A_u2f * A_dagger * input;
tau = A_u2tau * A_dagger * input;

% [force; tau] + ((A* (M\ A')) \ A)* inv(M) *input
%%
tic
%variable delta : 3n, t : 3
f = [zeros(3*n_am, 1); 1; 1;0];

A_eq = zeros(6, 3*n_am + 3);
for j = 1:n_am
    rj = [1; 0; 1];
    if j == 1
        rj = [1; 0; 2];
    end
    if j == 2
        rj = [2; 0; 3];
    end
    A_eq(1:6, 3*j-2:3*j) = [eye(3,3); S(rj)];
end
b_eq = [1,2,3,4,5,6]; % f_ext; re x f_ext + tau_ext

v = zeros(3*n_am, 1);
for j = 1:n_am
    v(3*j-2:3*j) = mass(j)* 9.81 * e_3;
end

% inf-norm of v - delta
A_ineq = [-1 * eye(3*n_am), -1 * ones(3*n_am, 1), zeros(3*n_am, 2) ...
          ;eye(3*n_am), -1 * ones(3*n_am, 1), zeros(3*n_am, 2)];

b_ineq = [-v; v];

% tilting angle constraint
A1 = zeros(n_am, 3*n_am);
A2 = zeros(n_am, 3*n_am);
tan_max = tan(1 / 180 * pi);
a1 = - e_1' + tan_max * e_3';
a2 = e_1' + tan_max * e_3';
for j = 1 : n_am
    Rj = eye(3);
    A1(j, 3*j-2:3*j) = a1 * Rj';
    A2(j, 3*j-2:3*j) = a2 * Rj';
end

A_ineq = [A_ineq; A1, zeros(n_am, 3);...
                  A2, zeros(n_am, 3)];

b_ineq = [b_ineq; A1 * v; A2 * v];

% inf-norm of internal force
% K_int = diag([1e1 * ones((n_am-1)*3, 1); 1e1 * ones((n_am-1)*3, 1)]);
K_int_inv = [1e1 \ ones((n_am-1)*3, 1); 1e1 \ ones((n_am-1)*3, 1)];
F_ext = [1;2; 3; zeros(3*n_am -3, 1); 0;0; 1; zeros(3*n_am -3, 1)]; %todo f_ext 0 ... ; t_ext
F_int0 =  A_dagger_2 * F_ext; %todo Adqd,  Cqd;
A_delta = A_dagger_2(:, 1:3*n_am);

A_ineq = [A_ineq; A_delta, zeros(6*(n_am-1), 1), -K_int_inv, zeros(6*(n_am-1), 1);...
                  -A_delta, zeros(6*(n_am-1), 1), -K_int_inv, zeros(6*(n_am-1), 1)];

b_ineq = [b_ineq; -F_int0; F_int0];

% inf-norm of difference
K_smooth_inv = 3e1 \ ones(n_am*3, 1);
delta_prev = zeros(3*n_am, 1);

A_ineq = [A_ineq; eye(3*n_am), zeros(3*n_am, 2), -K_smooth_inv;...
                  -eye(3*n_am), zeros(3*n_am, 2), -K_smooth_inv];

b_ineq = [b_ineq; delta_prev; -delta_prev];


lb = -inf * ones(3* n_am + 3, 1); 

ub = inf * ones(3* n_am + 3, 1);

options = optimoptions('linprog', ...
    'Algorithm', 'dual-simplex', ...  % 'interior-point', 'simplex', 'dual-simplex' 중 선택
    'Display', 'final', ...           % 'off', 'final', 'iter' 중 선택
    'MaxIterations', 100);             % 최대 반복 횟수 제한

[x, fval, exitflag, output] = linprog(f, A_ineq, b_ineq, A_eq, b_eq, lb, ub, options);
x(1:3*n_am)'
x(3*n_am+1 : end)'
toc