clear; close all;
addpath("../functions", "../../params")

%% Dynamics Parameters
params = define_params_ver2();
mu = params{3}; r = params{4};
thrust_limit = params{6};
B = [r r -r -r;  -r r r -r; mu -mu mu -mu; 1 1 1 1];

n_am = 4;
mass = 2.0 * ones(1,n_am); m_com = sum(mass);
k_spring = 300; k_s = [0, ones(1, n_am - 1) 0] * k_spring; % N/m , spring coeff
l1 = 0.4; l2 = 0.3;
pc_i = zeros(3, n_am);
for i = 1:n_am
    pc_i(1, i) = -(l1 + l2) * (i - (n_am+1)/2);
end
e_1 = [1;0;0]; e_2 = [0;1;0]; e_3 = [0;0;1];
g = 9.81;


a = 0.2; b = 0.2;c = 0.1;
J_quad = diag([b^2+c^2, a^2+c^2, a^2+b^2]) /5 * mass(1) * 0.7;

a = (l1+l2)/2; b = 0.05; c = 0.05;
J_tool = diag([b^2+c^2, a^2+c^2, a^2+b^2]) /5 * mass(1) * 0.3;


%% Pfaffian Constraints
M = zeros(n_am * 6, n_am * 6);
for j =1:n_am
    M(6*j-5: 6*j-3, 6*j-5: 6*j-3) = mass(j) * eye(3);
    M(6*j-2: 6*j, 6*j-2: 6*j) = J_quad + J_tool;
end

A = zeros((n_am-1)*6,n_am*6);
R1 = eye(3);
R2 = eye(3);

for j = 1:n_am - 1
    A(6*j-5: 6*j, 6*j-5: 6*j+6) = [eye(3), l2 * R1 *S(e_1), -eye(3), l1 * R2 *S(e_1); ...
        zeros(3,3), eye(3), zeros(3,3), -eye(3)];
end

A_dagger = A' * inv(A * inv(M) * A') * A * inv(M);

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

A_dagger * [ones(n_am*3, 1) * 1; ones(n_am*3,1)* 0];
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
        A_f_r(3*j-2:3*j, 3*j-5:3*j-3) = l1*S(Rj*e_1);
        A_tau_r(3*j-2:3*j, 3*j-5:3*j-3) = eye(3);
    else
        A_f_t(3*j-2:3*j, 3*j-5:3*j) = [eye(3), -eye(3)];
        A_f_r(3*j-2:3*j, 3*j-5:3*j) = [l1*S(Rj*e_1), -l2*S(Rj*e_1)];
        A_tau_r(3*j-2:3*j, 3*j-5:3*j) = [eye(3), -eye(3)];
    end
end
%
A_u2f = A_f_t(4:end, :)\A_dagger(4 : n_am*3, :);
A_u2tau = A_tau_r(4:end, :)\A_f_r(4:end, :)*A_u2f + A_tau_r(4:end, :)\ A_dagger((n_am+1)*3+1:end, :);
%%
A_u2f * A_dagger * [ones(n_am*3, 1) * 1; ones(n_am*3,1)* 1];