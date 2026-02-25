clear; close all;
addpath("../../params", "../../../casadi-3.6.7-windows64-matlab2018b")
import casadi.*
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

num_AMs = 4;
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

theta_min = -30 / 180 * pi * ones(num_AMs, 1);
theta_max = 30 / 180 * pi * ones(num_AMs, 1);

thrust_min = 1* -thrust_limit * ones(4*num_AMs, 1);
thrust_max = thrust_limit * ones(4*num_AMs, 1);
thrust_cen = (thrust_max + thrust_min)/2; % 4 x 1
thrust_tilde = (thrust_max(1) - thrust_min(1))/2; % scalar
%% fibonacci sphere
offset = m_c*  9.81 * [0; 0.0; 1.0];
num_points = 100; % Number of points
target_axis = offset;
cone_angle = 5; %degree
bound_ratio = 1.0;
U = sample_3d_cone_process(target_axis, cone_angle, num_points, bound_ratio)'; %  tau_y, fx, fz
beta = 1e-1; % 0.01 ~ 1.0
%%
N_k = num_points;    % 시나리오(k)의 개수
dim_f = 4*num_AMs;  % f 벡터의 차원
dim_u = 3;  % u 벡터의 차원
dim_th = num_AMs; % theta의 차원

f_min_val = -10 * ones(dim_f, 1);
f_max_val =  10 * ones(dim_f, 1);

% 예시 u_k 데이터 (무작위 생성)
u_data = U'; 

% ---------------------------------------------------------
% 2. 최적화 문제 정의 (Opti Stack)
% ---------------------------------------------------------
opti = casadi.Opti();

% --- High-Level 변수 ---
Theta = opti.variable(dim_th);

% B(theta) 행렬 정의 (예시: 회전 행렬 혹은 스케일링)
% 사용자의 문제에 맞는 B(theta) 함수로 대체하세요.
B = [];
for j = 1:num_AMs
    R = R_shape{j}' * [cos(Theta(j)), 0, -sin(Theta(j)); 0, 1, 0; sin(Theta(j)), 0, cos(Theta(j))]; % R_y(-theta)
    B_j = [R, hat(r_cj{j}) * R; zeros(3, 3), R] * [B_tau; zeros(2, 4); B_lambda];
    B = [B, B_j];
end
B_theta = [B(2, :); B(4, :); B(6, :)];

rho_all = []; % Objective 계산을 위해 rho들을 모을 리스트

% ---------------------------------------------------------
% 3. Loop over k (Low-Level 문제를 KKT로 변환하여 삽입)
% ---------------------------------------------------------
for k = 1:N_k
    u_k = u_data(:, k);
    
    % --- Low-Level Primal Variables ---
    rho = opti.variable(); 
    f = opti.variable(dim_f);
    
    % --- Low-Level Dual Variables ---
    nu = opti.variable(dim_u);      % Equality constraint (rho*u = B*f)
    lam_lb = opti.variable(dim_f);  % Lower bound inequality
    lam_ub = opti.variable(dim_f);  % Upper bound inequality
    
    % --- KKT Condition 1: Primal Feasibility ---
    % rho * u_k = B(theta) * f
    opti.subject_to( rho * u_k == mtimes(B_theta, f) );
    opti.subject_to( f_min_val <= f );
    opti.subject_to( f <= f_max_val );
    
    % --- KKT Condition 2: Stationarity (Lagrangian Gradient = 0) ---
    % Lagrangian L = -rho + nu'*(rho*u - B*f) + lam_ub'*(f - f_max) + lam_lb'*(f_min - f)
    
    % dL/drho = -1 + nu' * u_k = 0
    opti.subject_to( mtimes(nu', u_k) == 1 );
    
    % dL/df = -B(theta)' * nu + lam_ub - lam_lb = 0
    opti.subject_to( -mtimes(B_theta', nu) + lam_ub - lam_lb == 0 );
    
    % --- KKT Condition 3: Dual Feasibility ---
    opti.subject_to( lam_lb >= 0 );
    opti.subject_to( lam_ub >= 0 );
    
    % --- KKT Condition 4: Complementarity Slackness ---
    % 원래 조건: lam * slack == 0. 
    % 수치적 안정성을 위해 완화(Relaxation)된 조건을 주로 사용합니다 (NCP function 등).
    % 여기서는 간단히 직접 곱셈 형태로 구현하지만, 수렴이 어렵다면 epsilon을 사용하세요.
    
    % Strict complementarity (어려울 수 있음)
    % opti.subject_to( lam_lb .* (f - f_min_val) == 0 );
    % opti.subject_to( lam_ub .* (f_max_val - f) == 0 );
    
    % Relaxed complementarity (추천)
    epsilon = 1e-5;
    opti.subject_to( lam_lb .* (f - f_min_val) <= epsilon );
    opti.subject_to( lam_ub .* (f_max_val - f) <= epsilon );
    
    % 결과 수집
    beta = 0.1;
    rho_all = [rho_all; beta*rho];
end

% ---------------------------------------------------------
% 4. High-Level Objective
% ---------------------------------------------------------
% Maximize -logsumexp(rho_k)  <=> Minimize logsumexp(rho_k)
% CasADi에 logsumexp 내장 함수가 있습니다.
obj = logsumexp(rho_all); 

opti.minimize(obj);

% ---------------------------------------------------------
% 5. Solver 설정 및 실행
% ---------------------------------------------------------
% MPEC 문제는 비볼록(Non-convex)하므로 초기값이 중요합니다.
opti.set_initial(Theta, 0.5); 

p_opts = struct('expand', true); % 속도 향상
s_opts = struct('max_iter', 1000);
opti.solver('ipopt', p_opts, s_opts);

try
    sol = opti.solve();
    theta_opt = sol.value(Theta);
    rho_opt = sol.value(rho_all);
    fprintf('Optimal Theta: %f\n', theta_opt);
    disp('Optimal Rho values:');
    disp(rho_opt);
catch e
    disp('Solver failed or reached limits.');
    % 디버깅을 위해 현재 값 확인 가능
    theta_opt = opti.debug.value(Theta);
end