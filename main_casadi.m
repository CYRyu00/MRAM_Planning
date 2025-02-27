addpath("../casadi-3.6.7-windows64-matlab2018b")
import casadi.*

% Parameters
I = [16.571710 0.830806 0.718277;
     0.830806 16.655602 1.800197;
     0.718277 1.800197 29.261652]*10e-6;
I0 = I(2,2);
thrust_limit = 0.15 ;
m0 = 0.028;
mu = 0.005964552;
r = 0.092*sqrt(2)/4;
d = 0.1;
lg = 0.05; % custom
g=9.81;
m1 = 0.2; m2 = 0.1; lp = 0.2; l1 = lp+lg;
c_cart = 100e-2; 
c_pole = 100e-5; 

dt = 0.05;
%dt = 0.05;
nx = 4;
nu = 4;
N = 100;
%%
x = MX.sym('x');
dx = MX.sym('dx');
theta = MX.sym('theta');
dtheta = MX.sym('dtheta');
q = [x; theta; dx ;dtheta];

f1 = MX.sym('f1');
f2 = MX.sym('f2');
f3 = MX.sym('f3');
f4 = MX.sym('f4');
u = [f1;f2;f3;f4];
%%
    
% q1 q2 q1dot q2dot
M = [m1+m2+m0, (m2*lp+m0*l1)*cos(q(2));(m2*lp+m0*l1)*cos(q(2)) , m2*lp^2+m0*l1^2+I0];
C = [0 , -(m2*lp+m0*l1)*sin(q(2))*q(4); 0 0];
G = [0 ; (m2*lp+m0*l1)*g*sin(q(2))];

%Calculate Generalized force
A = [r r -r -r;-r r r -r;mu -mu mu -mu; zeros(2,4);ones(1,4)];
F_b = A*u;% [moment ; force]
R = [-sin(q(2)) 0 cos(q(2)); 0 -1 0; cos(q(2)) 0 sin(q(2))];
f_w = R*F_b(4:6);
J = [0 0;0 1;0 0; 1 l1*cos(q(2)); 0 0; 0 l1*sin(q(2))];
tau = J'*[F_b(1:3);f_w];


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
Q = diag([1,1,1,1])*0.1;
R = diag([1,1,1,1]);
x_init = MX.sym('x_init', nx);
X(:,1) = x_init;
x_0 = [0;0;0;0];
x_f = [3;pi/3;0;0];
u_max = thrust_limit*10e10;
u_min = -thrust_limit*0;

for k = 1:N
  % obj = obj + (X(:,k)-x_f)'*Q*(X(:,k)-x_f);
  Q = diag([0,0,1,1])*0.1;
  obj = obj + (X(:,k))'*Q*(X(:,k));

  obj = obj + U(:,k)'*R*U(:,k);
  g = [g; X(:,k+1) - f(X(:,k), U(:,k)); -X(:,k+1) + f(X(:,k), U(:,k))]; % dynamics constraint
  g = [g;  u_min - U(1,k);  u_min - U(2,k) ;  u_min - U(3,k);  u_min - U(4,k)];
  g = [g;  -u_max +  U(1,k);  -u_max +  U(2,k) ;  -u_max +  U(3,k);  -u_max + U(4,k)];
end
%g = [g; X(:,1) - x_0; - X(:,1) + x_0];
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

x0 = [0; 0; 0; 0];
u0 = [0; 0; 0; 0];

x_interp = zeros(N+1, 4);
vel_x = (x_f(1) - x_0(1))/(N*dt) ;
vel_theta = (x_f(2) - x_0(2))/(N*dt);

for i = 1:4
    x_interp(:, i) = linspace(x_0(i), x_f(i), N+1)';
end
for i = 1:(N+1)
    x_interp(i, 3:4) = [vel_x;vel_theta];

end

X0 = [reshape(x_interp',(N+1)*4,1);repmat(u0, N, 1)];
%[x_fmincon(1,:)';reshape(x_fmincon,N*4,1);reshape(u_fmincon,N*4,1)];
%[reshape(x_interp',(N+1)*4,1);repmat(u0, N, 1)];
%[repmat(x0, N+1, 1); repmat(u0, N, 1)]
%[x_fmincon(1,:)';reshape(x_fmincon,N*4,1);reshape(u_fmincon,N*4,1)];
%% solve
sol = solver('x0',  X0, ... 
             'p', x_0,...
             'lbx', -inf, 'ubx', inf,...
             'lbg', -inf, 'ubg', 0);


% Extract the optimal solution
solution = full(sol.x);
x_opt = reshape(solution(1:4*(N+1)), nx, N+1)';
u_opt = reshape(solution(4*(N+1)+1:end), nu, N)';
%% Plot the results
T = dt * N;
% state
t = (0:1:N)*dt;
figure(2)
subplot(2,2,1)
plot(t, x_opt(:,1),'.-')
ylabel("x1")
axis tight

subplot(2,2,2)
plot(t, x_opt(:,2),'.-')
ylabel("x2")

axis tight
subplot(2,2,3)
plot(t,x_opt(:,3),'.-')
ylabel("x3")
axis tight

subplot(2,2,4)
plot(t, x_opt(:,4),'.-')
ylabel("x4")
axis tight;

%% Inputs
figure(3)
t = (1:1:N)*dt;
subplot(2,2,1)
plot(t, u_opt(:,1),'.-')
ylabel("f1")
title("Inputs - AM1")
axis tight

subplot(2,2,2)
plot(t, u_opt(:,2),'.-')
ylabel("f2")

axis tight
subplot(2,2,3)
plot(t,u_opt(:,3),'.-')
ylabel("f3")

axis tight
subplot(2,2,4)
plot(t, u_opt(:,4),'.-')
ylabel("f4")
axis tight;
%%
taus = zeros(N,2);
Fs = zeros(N,6);

for i= 1:N
    tau_ = get_tau(x_opt(i,:)',u_opt(i,:)');
    taus(i,:)= tau_';
end
%%

figure(4);
subplot(2,1,1)
plot(t,taus(:,1),'.-');

title("tau1")
subplot(2,1,2)
plot(t,taus(:,2),'.-');

title("tau2")