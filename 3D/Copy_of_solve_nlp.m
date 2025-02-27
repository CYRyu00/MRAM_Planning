function [x_opt, u_opt, optimal_value,exit_flag,processing_time] ....
          =  solve_Nlp(params,shape,num_AMs,dh,gravity, mass,inertia,r_i_ci,qo_desired,tau_min, tau_max, u_min,u_max,x_0,x_f,X_init_guess,dt,N)
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
%F_ext = [0;0;0;0;0;0]; tau = [0;0;0;0];
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
Q1 = diag([0,0,0,0,1,1,1,1])*0.1; Q2 =diag([1,1]);
R = diag(ones(nu,1))*1;


for k = 1:N
  % obj = obj + (X(:,k)-x_f)'*Q*(X(:,k)-x_f);
  obj = obj + X(:,k)'*Q1*X(:,k) + (X(1:2,k) - qo_desired(:,k))'*Q2*(X(1:2,k) - qo_desired(:,k));
  obj = obj + U(:,k)'*R*U(:,k);

  % dynamics equality constraint
  g = [g; X(:,k+1) - f(X(:,k), U(:,k))]; 
  
end
%g = [g; X(:,N+1) - x_f; - X(:,N+1) + x_f];
g = [g; X(:,1) - x_0];
g = [g; X(:,N+1) - x_f];

%%
nlp_prob = struct('x', opt_variables, 'f', obj, 'g', g);
nlp_opts = struct;
nlp_opts.ipopt.print_level = 1;
nlp_opts.ipopt.tol = 1e-6;
nlp_opts.ipopt.max_iter = 300;
nlp_opts.ipopt.mu_strategy = 'adaptive';
nlp_opts.ipopt.linear_solver = 'mumps';
%nlp_opts.ipopt.jacobian_approximation = 'exact'; 
nlp_opts.ipopt.hessian_approximation = 'limited-memory';
nlp_opts.print_time = 1;

solver = nlpsol('solver', 'ipopt', nlp_prob, nlp_opts);
%% solve
LBx = [ repmat([ -0.1*pi; -pi/3; -pi/2 ; -pi/2; -inf; -inf; -inf; -inf], N+1, 1); ...
        repmat([tau_min; tau_min; u_min*ones(nu-2,1)], N, 1)];
UBx = [ repmat([ pi/2; pi/3; pi/2 ; pi/2; inf; inf; inf; inf], N+1, 1); ...
        repmat([tau_max; tau_max; u_max*ones(nu-2,1)], N, 1)];

sol = solver('x0',  X_init_guess, ... 
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


n = 4;
% dynamic parameters of each module 
params = define_params();
m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5}; thrust_limit= params{6};

dh = [0,0,0.95,0;   % [alpha, a, d, theta]
      -pi/2, 0.9 , 0,0;
      0,-0.1,0.23,pi/2;
      pi/2,0,0,-pi/2;
      pi/2,0,0,0];
gravity = [0;0;-9.81];

m=7;
%all_shapes = generate_all_shapes(m);

% Show example rigidbodytree
shape = all_shapes{10}{10}; num_AMs= 10;
[AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);
mass =  {10, 1, 0.5, AM_mass};
inertia = {eye(3)*1, eye(3)*0.1, eye(3)*0.1, AM_inertia, zeros(3,3)};
r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};

do_view = 0; q=[pi/6;-pi/3;pi/3;-pi/4];
%robot = generate_door_ver2(n,dh,r_i_ci, gravity, mass,inertia, do_view,q);

% NLP parameters
x_0 = [0;0;0;0;0;0;0;0];
x_f = [pi/4;-pi/6;pi/6;0;0;0;0;0];
thrust_scale = 3;
u_max = thrust_limit *thrust_scale;
u_min = thrust_limit *(-thrust_scale);
tau_scale = 3;
tau_min = -0.2 *tau_scale ; 
tau_max =  0.2 *tau_scale ;
dt = 0.1;
N = 100;

qo_desired = zeros(2,N+1);
T = N*dt; t0 = 1; %1sec
for j=1:N+1
    t = (j-1)*dt;
    if t < t0 
        qo_desired(:,j) = [0 ; pi/6*(cos(pi/t0*t) - 1)/2];
    elseif t < T-t0
        qo_desired(:,j) = [(x_f(1)-x_0(1))/(T-t0)*(t-t0) ;-pi/6];
    else
        qo_desired(:,j) = [(x_f(1)-x_0(1))/(T-t0)*(t-t0) ;pi/6*(cos(pi/t0*t -pi/t0*T) - 1)/2];
    end
end


x_interp = zeros(N+1, 8);
vel_x1 = (x_f(1) - x_0(1))/(N*dt); vel_x2 = (x_f(2) - x_0(2))/(N*dt);
vel_x3 = (x_f(3) - x_0(3))/(N*dt);vel_x4 = (x_f(4) - x_0(4))/(N*dt);
for j = 1:8
    x_interp(:, j) = linspace(x_0(j), x_f(j), N+1)';
end
for j = 1:(N+1)
    x_interp(j, 5:8) = [vel_x1;vel_x2;vel_x3;vel_x4];
end


nu = 2 + num_AMs*4; zero_us = zeros(nu,1); 
X_init_guess = [reshape(x_interp',(N+1)*8,1);repmat(zero_us, N, 1)];

[ x_opt, u_opt, optimal_value,exit_flag,processing_time] ....
          = solve_Nlp(params,shape,num_AMs,dh, gravity,mass,inertia,r_i_ci, ...
                      qo_desired, tau_min, tau_max,u_min,u_max,x_0,x_f,X_init_guess,dt,N);
fprintf("Exit Flag: %d\n" ,exit_flag)
%%
close all
figure;
time = (1:1:N+1)*dt;
subplot(3,1,1)
plot(time, x_opt(:,1:4))
hold on
plot(time, qo_desired, "--")
legend("q1","q2","q3","q4","q1\_desired","q2\_desired")
title("states")
axis tight;

subplot(3,1,2)
plot(time(1:end-1),u_opt(:,1:2))
legend("u1","u2")
title("motor inputs")
axis tight;

wrench = zeros(N,6);
for i=1:N
    wrench(i,:) = map_u2wrench_double( u_opt(i,3:end)', shape , mu , r , d);
end
subplot(3,1,3)
plot(time(1:end-1),wrench)
legend("m_x","m_y","m_z","f_x","f_y","f_z")
title("Wrench exp. AM frame")
axis tight;

%%
do_view = 0 ;
robot = generate_door_ver2(n,dh,r_i_ci, gravity, mass,inertia, do_view,q);
slow_factor = 1;
force_scale = 0.1;
plot_tree(robot, dh, params, x_opt,u_opt, dt,N,slow_factor, force_scale, shape)
