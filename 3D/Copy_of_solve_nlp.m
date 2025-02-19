function [x_opt, u_opt, optimal_value,exit_flag,processing_time] ....
          =  solve_Nlp(params,shape,num_AMs,dh,gravity, mass,inertia,r_i_ci,tau_min, tau_max, u_min,u_max,x_0,x_f,X_init_guess,dt,N)
addpath("../../casadi-3.6.7-windows64-matlab2018b")
import casadi.*

n=4;
m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5};

% 
nx = n*2;
nu = 2 + num_AMs*4;

x = MX.sym('x', nx, 1); % [q; qd];
u = MX.sym('u', nu, 1);
%% Dynamics
F_ext = [0;0;0;0;0;0]; tau = [0;0;0;0];
tau = [0;0;u(1);u(2)];
F_ext = map_u2wrench( u(3:end), shape , mu , r , d);

qdd = FD_ver2(n, dh, mass, inertia, r_i_ci, gravity, x(1:4), x(5:8), tau, F_ext);

x_dot = [x(5:8);qdd];
%x_damp = [0;0;-c_cart*q(3); -c_pole*q(4)];
x_damp = [0;0;0;0;0;0;0;0];
x_next = x + dt * ( x_dot + x_damp) ; 
% Create CasADi function for the dynamics
f = Function('f', {x, u}, {x_next});
%%
U = MX.sym('U', nu, N);
X = MX.sym('X', nx, N+1);
opt_variables = [reshape(X, nx*(N+1), 1); reshape(U, nu*N, 1)];
%%
obj = 0;
g = [];
Q = diag(ones(nx,1))*0.1;
R = diag(ones(nu,1))*1;
x_init = MX.sym('x_init', nx);
X(:,1) = x_init;
%u_max = thrust_limit *1;
%u_min = thrust_limit *(-0);

for k = 1:N
  % obj = obj + (X(:,k)-x_f)'*Q*(X(:,k)-x_f);
  Q = diag([0,0,0,0,1,1,1,1])*0.1;
  obj = obj + X(:,k)'*Q*X(:,k);
  obj = obj + U(:,k)'*R*U(:,k);

  % dynamics constraint
  %g = [g; X(:,k+1) - f(X(:,k), U(:,k)); -X(:,k+1) + f(X(:,k), U(:,k))]; 
  g = [g; X(:,k+1) - f(X(:,k), U(:,k))]; 
  %Inequality Constraints
  %g = [g;  u_min*ones(nu-2,1) - U(3:end,k)];
  %g = [g; -u_max*ones(nu-2,1) + U(3:end,k)];
  
end
%g = [g; X(:,N+1) - x_f; - X(:,N+1) + x_f];
g = [g; X(:,N+1) - x_f];

%%
nlp_prob = struct('x', opt_variables, 'f', obj, 'g', g, 'p', x_init);
nlp_opts = struct;
nlp_opts.ipopt.print_level = 1;
nlp_opts.ipopt.tol = 1e-6;
nlp_opts.ipopt.mu_strategy = 'adaptive';
nlp_opts.ipopt.linear_solver = 'mumps';
%nlp_opts.ipopt.jacobian_approximation = 'exact'; 
nlp_opts.ipopt.hessian_approximation = 'limited-memory';
nlp_opts.print_time = 1;

solver = nlpsol('solver', 'ipopt', nlp_prob, nlp_opts);
%% solve


LBx = [ repmat([ -0.1*pi; -pi/2; -pi ; -pi; -inf; -inf; -inf; -inf], N+1, 1); ...
        repmat([tau_min; tau_min; u_min*ones(nu-2,1)], N, 1)];
UBx = [ repmat([ pi/2; pi/2; pi ; pi; inf; inf; inf; inf], N+1, 1); ...
        repmat([tau_max; tau_max; u_max*ones(nu-2,1)], N, 1)];

sol = solver('x0',  X_init_guess, ... 
             'p', x_0,...
             'lbx', LBx, 'ubx', UBx,...
             'lbg', 0, 'ubg', 0);


% Extract the optimal solution
solution = full(sol.x);
x_opt = reshape(solution(1:nx*(N+1)), nx, N+1)';
u_opt = reshape(solution(nx*(N+1)+1:end), nu, N)';
% Get the informations
optimal_value = full(sol.f);
exit_flag = solver.stats.success;
processing_time = solver.stats.t_wall_total;
end

n=4;
dh = [0,0,0.95,0;   % [alpha, a, d, theta]
      -pi/2, 0.9 , 0,0;
      0,-0.1,0.23,pi/2;
      pi/2,0,0,-pi/2;
      pi/2,0,0,0];
gravity = [0;0;-9.81];

params = define_params();
m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5};thrust_limit= params{6};

m=6; num_AMs = m;
all_shapes = generate_all_shapes(m);
shape = all_shapes{6}{10};
[AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);
mass = {10, 1, 0.5, AM_mass};
inertia = {eye(3)*1, eye(3)*2, eye(3)*3, AM_inertia,};
r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)]};

x_0 = [0;0;0;0;0;0;0;0];
x_f = [pi/3;-pi/6;pi/6;0;0;0;0;0];
thrust_scale = 3;
u_max = thrust_limit *thrust_scale;
u_min = thrust_limit *(-thrust_scale);
tau_min = -0.75; tau_max = 0.75;
dt = 0.1;
N = 100;

x_interp = zeros(N+1, 8);
vel_x1 = (x_f(1) - x_0(1))/(N*dt); vel_x2 = (x_f(2) - x_0(2))/(N*dt);
vel_x3 = (x_f(3) - x_0(23))/(N*dt);vel_x4 = (x_f(4) - x_0(4))/(N*dt);
for j = 1:8
    x_interp(:, j) = linspace(x_0(j), x_f(j), N+1)';
end
for j = 1:(N+1)
    x_interp(j, 5:8) = [vel_x1;vel_x2;vel_x3;vel_x4];
end

nu = 2 + num_AMs*4; zero_us = zeros(nu,1); 
X_init_guess = [reshape(x_interp',(N+1)*8,1);repmat(zero_us, N, 1)];


[x_opt, u_opt, optimal_value,exit_flag,processing_time] ....
          = solve_Nlp(params,shape,num_AMs,dh, gravity,mass,inertia,r_i_ci, ...
                      tau_min, tau_max,u_min,u_max,x_0,x_f,X_init_guess,dt,N);
%%
close all
time = (1:1:N+1)*dt;
subplot(2,1,1)
plot(time, x_opt(:,1:4))
legend("q1","q2","q3","q4")
axis tight;
subplot(2,1,2)
plot(time(1:end-1),u_opt(:,1:6))
legend("u1","u2","u3","u4","u5","u6")
axis tight;