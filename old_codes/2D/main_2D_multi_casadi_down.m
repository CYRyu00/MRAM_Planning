addpath("../../casadi-3.6.7-windows64-matlab2018b")
import casadi.*

% Parameters
I = [16.571710 0.830806 0.718277;
     0.830806 16.655602 1.800197;
     0.718277 1.800197 29.261652]*10e-6;
I0 = I(2,2);
thrust_limit = 0.15 ;
m0 = 0.028;
mu = 0.005964552;
r = 0.092*sqrt(2)/4;%0.032
d = 0.1;
lg = 0.08; % custom
g=9.81;
m1 = 0.2; m2 = 0.1; lp = 0.2; l1 = lp+lg;
c_cart = 100e-2; 
c_pole = 100e-5; 

global params
params = [m1, m2, lp, lg, m0, I(2,2),mu,r];
dt = 0.05;
N = 100;

%distance from box to CoM of each AM

Max_serial_num = 4;
L = (lp+lg):d:(lp+lg + (Max_serial_num-1)*d);

num_up = [3 3 1 0];  %[3 3 1];
num_down = [0 0 0 0]; %[0 2 0];
num_AMs = sum(num_up)+sum(num_down);
nx = 4;
nu = ( length(num_up) + length(num_down) )*2;


x = MX.sym('x');
dx = MX.sym('dx');
theta = MX.sym('theta');
dtheta = MX.sym('dtheta');
q = [x; theta; dx ;dtheta];

u = MX.sym('u', ( length(num_up) + length(num_down) )*2, 1);
%% Dynamics
% q1 q2 q1dot q2dot

term1 = m1 + m2;
term2 = m2*lp;
term3 = m2*lp^2;
for i=1:length(num_up)
    term1 = term1 + m0 *(num_up(i) + num_down(i));
    term2 = term2 + m0*L(i) *(num_up(i) + num_down(i));
    term3 = term3 + (m0*L(i)^2 + I0) *(num_up(i) + num_down(i));
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
for i=1:length(num_down)
    F_b = A*u(length(num_up)*2 + 2*i-1: length(num_up)*2 +2*i)*num_down(i);% [moment ; force]
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
R = diag([reshape([num_up; num_up], 1, []),reshape([num_down; num_down], 1, [])])*2;
x_init = MX.sym('x_init', nx);
X(:,1) = x_init;
x_0 = [0;0;0;0];
x_f = [3;pi/3;0;0];
u_max = thrust_limit *10e10;
u_min = thrust_limit *(-0);

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
  for i=1:length(num_down)
      g = [g;  u_min*ones(2,1) - U(length(num_up)*2+ 2*i-1:length(num_up)*2+ 2*i,k)];
      g = [g; -u_max*ones(2,1) + U(length(num_up)*2+ 2*i-1:length(num_up)*2+ 2*i,k)];
  end
   
end
g = [g; X(:,N+1) - x_f; - X(:,N+1) + x_f];

%%
nlp_prob = struct('x', opt_variables, 'f', obj, 'g', g, 'p', x_init);
nlp_opts = struct;
nlp_opts.ipopt.print_level = 3;
nlp_opts.ipopt.tol = 1e-6;
nlp_opts.ipopt.mu_strategy = 'monotone';
nlp_opts.ipopt.linear_solver = 'mumps';
nlp_opts.print_time = 1;

solver = nlpsol('solver', 'ipopt', nlp_prob, nlp_opts);

zero_xs = zeros(length(x),1);
zero_us = zeros(length(u),1);

x_interp = zeros(N+1, 4);
vel_x = (x_f(1) - x_0(1))/(N*dt) ;
vel_theta = (x_f(2) - x_0(2))/(N*dt);

for i = 1:4
    x_interp(:, i) = linspace(x_0(i), x_f(i), N+1)';
end
for i = 1:(N+1)
    x_interp(i, 3:4) = [vel_x;vel_theta];

end

X_init_guess = [reshape(x_interp',(N+1)*4,1);repmat(zero_us, N, 1)];
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
%% Plot
close all
plot_results(x_opt, u_opt,dt,N);
%% Video
slow_factor = 2;
plot_tree(x_opt,u_opt, nx,nu,dt,N,slow_factor)