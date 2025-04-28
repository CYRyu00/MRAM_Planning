function [x_opt, u_opt, optimal_value,exit_flag,processing_time] ....
          = solve_nlp(num_up,L,u_min,u_max,x_0,x_f,X_init_guess,dt,N)
addpath("../../casadi-3.6.7-windows64-matlab2018b")
import casadi.*

global params
m1 = params(1);
m2 = params(2);
lp = params(3);
lg = params(4);
m0 = params(5);
I0 = params(6);
mu = params(7);
r = params(8);
d = params(9);
g = params(10);
c_cart = params(11);
c_pole = params(12);


%dt = 0.05;
%N = 100;

%distance from box to CoM of each AM

%Max_serial_num = 4;
%L = (lp+lg):d:(lp+lg + (Max_serial_num-1)*d);


%num_up = [3 1 4];  %[3 3 1];

num_AMs = sum(num_up);
nx = 4;
nu = length(num_up)*2;


x = MX.sym('x');
dx = MX.sym('dx');
theta = MX.sym('theta');
dtheta = MX.sym('dtheta');
q = [x; theta; dx ;dtheta];

u = MX.sym('u', length(num_up)*2, 1);
%% Dynamics
% q1 q2 q1dot q2dot

term1 = m1 + m2;
term2 = m2*lp;
term3 = m2*lp^2;
for i=1:length(num_up)
    term1 = term1 + m0 *num_up(i) ;
    term2 = term2 + m0*L(i) *num_up(i);
    term3 = term3 + (m0*L(i)^2 + I0) *num_up(i);
end
% q1 q2 q1dot q2dot
M = [term1, term2*cos(q(2));term2*cos(q(2)) , term3];
C = [0 , -term2*sin(q(2))*q(4); 0 0];
G = [0 ; term2*g*sin(q(2))];
    

%Calculate Generalized force
A = [0 0;-2*r 2*r; 0 0; 0 0; 0 0 ; 2 2];
tau = 0;

for i=1:length(num_up)
    F_b = A*u(2*i-1:2*i)*num_up(i);% [moment ; force]
    R = [-sin(q(2)) 0 cos(q(2)); 0 -1 0; cos(q(2)) 0 sin(q(2))];
    f_w = R*F_b(4:6);
    J = [0 0;0 1;0 0; 1 L(i)*cos(q(2)); 0 0; 0 L(i)*sin(q(2))];
    tau = tau + J'*[F_b(1:3);f_w];
end

q_dot = [q(3:4); M\(tau - C*q(3:4)-G)];
q_damp = [0;0;-c_cart*q(3); -c_pole*q(4)];
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
Q = diag(ones(length(x),1))*0.1;
R = diag(reshape([num_up; num_up], 1, []))*2;
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
solution = full(sol.x);
x_opt = reshape(solution(1:4*(N+1)), nx, N+1)';
u_opt = reshape(solution(4*(N+1)+1:end), nu, N)';
% Get the informations
optimal_value = full(sol.f);
exit_flag = solver.stats.success;
processing_time = solver.stats.t_wall_total;
