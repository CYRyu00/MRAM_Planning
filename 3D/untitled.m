%%
setenv('CASADIPATH', 'C:\SNU\졸업논문\Matlab_project\casadi-windows-matlabR2016a-v3.5.5');
setenv('PATH', [getenv('PATH') ';C:\SNU\졸업논문\Matlab_project\casadi-windows-matlabR2016a-v3.5.5']);
disp(getenv('CASADIPATH'))  % CASADIPATH 확인
disp(getenv('PATH'))        % PATH 확인

%%
addpath("../../casadi-windows-matlabR2016a-v3.5.5")
import casadi.*
disp('CasADi 버전:');
disp(casadi.CasadiMeta.version());
methods(casadi.MX)
%%
A = MX.ones(3,1);
B = MX.ones(3,1);
cross(A,B)
sin(A)
%%
dll_path = fullfile('C:\SNU\졸업논문\Matlab_project\casadi-windows-matlabR2016a-v3.5.5', 'libcasadi_expm_slicot.dll');

try
    loadlibrary(dll_path, '');
    disp('✅ DLL 파일 로드 성공!');
catch ME
    disp('❌ DLL 로드 실패:');
    disp(ME.message);
end

%%
A_test = MX.eye(3); % 3x3 단위 행렬
try
    expm(A_test) % CasADi의 expm 실행 테스트
    disp('CasADi의 expm 함수가 정상적으로 작동합니다.')
catch ME
    disp('expm 실행 오류:')
    disp(ME.message)
end
%%
import casadi.*

m=8;num_AMs = m;
all_shapes = generate_all_shapes(m);
shape = all_shapes{8}{47};

params = define_params();
m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5};
[AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);

n=4;
M = cell(1,n+1);
A = cell(1,n);
w = { [0;0;1], [0;1;0], [1;0;0], [0;1;0] };
r = { [-0.5;-0.05;-1], [0.05;-0.07;0.0], [0.03;0;0], -AM_com };
for i=1:n
    A{i} = [ w{i}; -hat(w{i})*r{i}];
end

%M{i} is M_i,i-1
M{1} = inv( [ eye(3,3), [0.5;0.05;1] ; 0,0,0,1 ]); % base frame to door 
M{2} = inv( [ eye(3,3), [0.4;0.1;-0.05] ; 0,0,0,1 ]); % door to handle&gripper
M{3} = inv( [ [0,1,0; -1,0,0; 0,0,1] , [-0.03; 0.1 ; 0.0] ; 0,0,0,1 ]); % handle&gripper to link
M{4} = inv( [ eye(3,3), AM_com+[-0.05;0;0] ; 0,0,0,1 ]);% link to AM
M{5} = inv( [ eye(3,3), -AM_com ; 0,0,0,1 ]);% AM to core module

mass = {20, 3, 0.5, AM_mass};
inertia = {eye(3)*1, eye(3)*2, eye(3)*3, AM_inertia};
g = [0;0;-9.81];

% 
nx = n*2;
nu = num_AMs*4;
q1 = MX.sym('q1');
q2 = MX.sym('q2');
q3 = MX.sym('q3');
q4 = MX.sym('q4');
q1d = MX.sym('q1d');
q2d = MX.sym('q2d');
q3d = MX.sym('q3d');
q4d = MX.sym('q4d');
q = [q1;q2;q3;q4;q1d;q2d;q3d;q4d];

u = MX.sym('u', num_AMs*4, 1);
%% Dynamics
F_tip = [0;0;0;0;0;0]; tau = [0;0;0;0];
expm( -hat6(A{i}) *q(1))
q_dot = FD(n, A, M, mass, inertia, F_tip, q(1:4), q(5:8), tau, g);
%q_damp = [0;0;-c_cart*q(3); -c_pole*q(4)];
q_next = q + dt * (q_dot + q_damp) ; 
% Create CasADi function for the dynamics
f = Function('f', {q, u}, {q_next});
%%
U = MX.sym('U', nu, N);
X = MX.sym('X', nx, N+1);
opt_variables = [reshape(X, nx*(N+1), 1); reshape(U, nu*N, 1)];
%%
obj = 0;
g = [];
Q = diag(ones(length(X(:,k)),1))*0.1;
R = diag(ones(length(U(:,k)),1))*2;
x_init = MX.sym('x_init', nx);
X(:,1) = x_init;
%x_0 = [0;0;0;0];
%x_f = [3;pi/3;0;0];
%u_max = thrust_limit *1;
%u_min = thrust_limit *(-0);

for k = 1:N
  % obj = obj + (X(:,k)-x_f)'*Q*(X(:,k)-x_f);
  Q = diag([0,0,1,1])*0.1;
  obj = obj + (X(:,k))'*Q*(X(:,k));
  obj = obj + U(:,k)'*R*U(:,k);

  % dynamics constraint
  g = [g; X(:,k+1) - f(X(:,k), U(:,k)); -X(:,k+1) + f(X(:,k), U(:,k))]; 
  %Inequality Constraints
  for i=1:length(num_up)
      g = [g;  u_min*ones(2,1) - U(2*i-1:2*i,k)];
      g = [g; -u_max*ones(2,1) + U(2*i-1:2*i,k)];
  end
end
g = [g; X(:,N+1) - x_f; - X(:,N+1) + x_f];

%%
nlp_prob = struct('x', opt_variables, 'f', obj, 'g', g, 'p', x_init);
nlp_opts = struct;
nlp_opts.ipopt.print_level = 1;
nlp_opts.ipopt.tol = 1e-6;
nlp_opts.ipopt.mu_strategy = 'monotone';
nlp_opts.ipopt.linear_solver = 'mumps';
nlp_opts.print_time = 1;

solver = nlpsol('solver', 'ipopt', nlp_prob, nlp_opts);

x_interp = zeros(N+1, 4);
vel_x = (x_f(1) - x_0(1))/(N*dt) ;
vel_theta = (x_f(2) - x_0(2))/(N*dt);

for i = 1:4
    x_interp(:, i) = linspace(x_0(i), x_f(i), N+1)';
end
for i = 1:(N+1)
    x_interp(i, 3:4) = [vel_x;vel_theta];

end

%X_init_guess = [reshape(x_interp',(N+1)*4,1);repmat(zero_us, N, 1)];
%[x_fmincon(1,:)';reshape(x_fmincon,N*4,1);reshape(u_fmincon,N*4,1)];
%[reshape(x_interp',(N+1)*4,1);repmat(u0, N, 1)];
%[repmat(zero_xs, N+1, 1); repmat(zero_us, N, 1)]
%[x_fmincon(1,:)';reshape(x_fmincon,N*4,1);reshape(u_fmincon,N*4,1)];
%% solve
sol = solver('x0',  X_init_guess, ... 
             'p', x_0,...
             'lbx', -inf, 'ubx', inf,...
             'lbg', -inf, 'ubg', 0);


% Extract the optimal solution
solution = full(sol.x);a
x_opt = reshape(solution(1:4*(N+1)), nx, N+1)';
u_opt = reshape(solution(4*(N+1)+1:end), nu, N)';
% Get the informations
optimal_value = full(sol.f);
exit_flag = solver.stats.success;
processing_time = solver.stats.t_wall_total;
