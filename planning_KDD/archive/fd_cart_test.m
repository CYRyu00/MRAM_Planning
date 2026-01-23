clear; close all;
addpath("../params")
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

dt = 0.1; N = 10/dt;

mass_obj = 10;
k_obj = 20;
b_obj = 30; % crt:28.28
n_o = 3;

n_am = 4;
mass_am = 2.0 * ones(1,n_am); m_com = sum(mass_am);
l1 = 0.4; l2 = 0.3;
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

phi_1 = 0 / 180 * pi;
p_1 = r_g -  l1 * Ry(phi_1) * e_1 + [1; 0; 0];
pd_1 = [0;0;0];
w_1 = [0; 1; 0]*0;

qd = zeros(n_o + n_am*6, 1);
p_o = p_1 - r_g + l1 * Ry(phi_1) * e_1;
qd(1:3) = pd_1 + l1 * Ry(phi_1) * S(w_1) * e_1;

phi(1) = phi_1;
qd(4:6) = pd_1;
pd_prev = pd_1;
qd(n_o + (n_am + 1)*3 - 2 : n_o + (n_am + 1)*3) = w_1;

for j = 2:n_am
    phi(j) = phi_1 + shape_pitch(j);
    qd(n_o + (n_am + j)*3 - 2 : n_o + (n_am + j)*3) = w_1;
    qd(1 + 3*j : 3 + 3*j) = pd_prev - l2*Ry(phi(j-1))*S(w_1)*e_1 - l1*Ry(phi(j))*S(w_1)*e_1 ;
end

U = [randn((n_am)*3, 1) * 3; randn((n_am)*3, 1) * 1];
U = [repmat([10;1;1], n_am, 1); randn((n_am)*3, 1) * 0];

[qdd, f_in, tau_in] = FD_cart(p_o, phi, qd, U, M, k_obj, b_obj, l1, l2, n_o, n_am);
[qdd2, f_in2, tau_in2] = FD_cart_ver2(p_1, phi_1, pd_1, w_1, shape_pitch, r_g, U, M, k_obj, b_obj, l1, l2, n_o, n_am);
qdd2
%%
function out = Ry(q)
    out = [cos(q), 0, sin(q);
           0, 1, 0;
           -sin(q), 0, cos(q)];
end