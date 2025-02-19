function [x_opt, u_opt, optimal_value,exit_flag,processing_time] ....
          =  solve_nlp(params,shape,num_AMs,dh,gravity, mass,inertia,r_i_ci,tau_min, tau_max, u_min,u_max,x_0,x_f,X_init_guess,dt,N)
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
