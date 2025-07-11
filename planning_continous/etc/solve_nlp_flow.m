function [rho_opt, term_opt,flow_opt1,flow_opt2,flow_opt3,flow_opt4, x_opt, u_opt, optimal_value,exit_flag,processing_time] ....
          =  solve_NLP(params,num_AMs,K,L,core,dh,gravity,qo_desired,tau_min, tau_max, u_min,u_max,x_0,x_f,X_init_guess,dt,N)
addpath("../../casadi-3.6.7-windows64-matlab2018b")
import casadi.*

n=4;
m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5}; kt= params{7};c_1=params{8};c_2=params{9};

% 

nx = n*2;
nu = 2 + K*L*4;

x = MX.sym('x', nx, 1); % [q; qd];
u = MX.sym('u', nu, 1);
rho = MX.sym('rho', K,L); %shape

flow = cell(4,1);
for d = 1:4
    flow{d} = MX.sym(sprintf('flow_%d',d), K, L);
end
term = MX.sym('term', K, L);

%% Dynamics
%F_ext = [0;0;0;0;0;0]; tau = [0;0;0;0];

[AM_com, AM_mass, AM_inertia] = get_inertia(rho,K,L, core ,m0, I0, d);
mass =  {10, 1, 0.5, AM_mass};
inertia = {eye(3)*1, eye(3)*0.01, eye(3)*0.01, AM_inertia, zeros(3,3)};
r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};

tau = [-c_1*x(5);(-c_2*x(6) -kt*x(2) + mass{2}*0.80*9.81*0.040 ); u(1); u(2)] ;
F_ext = map_u2wrench( u, rho,K,L, core , mu , r , d);

qdd = FD_ver2(n, dh, mass, inertia, r_i_ci, gravity, x(1:4), x(5:8), tau, F_ext);

x_dot = [x(5:8);qdd];
x_next = x + dt *  x_dot ; 
% Create CasADi function for the dynamics
f = Function('f', {x, u , rho}, {x_next});
% Flow per direction: 1=Up, 2=Down, 3=Left, 4=Right

%%
U = MX.sym('U', nu, N);
X = MX.sym('X', nx, N+1);
opt_variables = [reshape(rho, K*L,1); reshape(flow{1},K*L,1); reshape(flow{2},K*L,1); reshape(flow{3},K*L,1); reshape(flow{4},K*L,1); reshape(term,K*L,1);...
                reshape(X, nx*(N+1), 1); reshape(U, nu*N, 1)];
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
  g = [g; X(:,k+1) - f(X(:,k), U(:,k),rho) ; -X(:,k+1) + f(X(:,k), U(:,k),rho)]; 
end

g = [g; rho(core(1),core(2))-1; -rho(core(1),core(2))+1];

eps = 0.001;

for i_=1:K
    for j_=1:L

        flow_in = MX(0); flow_out = MX(0);

        % In & out flows
        if i_ < K
            flow_in = flow_in + flow{1}(i_+1,j_); % From below
            flow_out = flow_out + flow{2}(i_,j_); % Down
            g = [g; -flow{1}(i_+1,j_)*flow{2}(i_,j_)-eps;flow{1}(i_+1,j_)*flow{2}(i_,j_)-eps];
        end  
        if i_ > 1
            flow_in = flow_in + flow{2}(i_-1,j_); % From above
            flow_out = flow_out + flow{1}(i_,j_); % Up
            g = [g; -flow{2}(i_-1,j_)*flow{1}(i_,j_)-eps;flow{2}(i_-1,j_)*flow{1}(i_,j_)-eps];
        end  
        if j_ < L
            flow_in = flow_in + flow{3}(i_,j_+1); % From right
            flow_out = flow_out + flow{4}(i_,j_); % Right
            g = [g; -flow{3}(i_,j_+1)*flow{4}(i_,j_)-eps;flow{3}(i_,j_+1)*flow{4}(i_,j_)-eps];
        end  
        if j_ > 1
            flow_in = flow_in + flow{4}(i_,j_-1); % From left
            flow_out = flow_out + flow{3}(i_,j_); % Left
            g = [g; -flow{4}(i_,j_-1)*flow{3}(i_,j_)-eps;flow{4}(i_,j_-1)*flow{3}(i_,j_)-eps];
        end  

        if i_ == 1, g=[g; flow{1}(i_,j_);-flow{1}(i_,j_)]; end  % Up
        if i_ == K, g=[g; flow{2}(i_,j_);-flow{2}(i_,j_)]; end  % Down
        if j_ == 1, g=[g; flow{3}(i_,j_);-flow{3}(i_,j_)]; end  % Left
        if j_ == L, g=[g; flow{4}(i_,j_);-flow{4}(i_,j_)]; end  % Right
        
        if i_==core(1)&j_==core(2)
            flow_in = MX(1);
            g = [g; term(i_,j_);-term(i_,j_) ];
        else
            g = [g; rho(i_,j_)-fmin(1,flow_in);-rho(i_,j_)+fmin(1,flow_in)];
        end
        %flow_in = flow_out;
        g = [g; (flow_in*(1-term(i_,j_)) - flow_out);(-flow_in*(1-term(i_,j_)) + flow_out)];
        
        for d = 1:4
            g = [g; flow{d}(i_,j_)*(flow{d}(i_,j_)-1)-eps; -flow{d}(i_,j_)*(flow{d}(i_,j_)-1)-eps];
            g = [g; -flow{d}(i_,j_); flow{d}(i_,j_)-1];
        end
        %g = [g; rho(i_,j_)*(rho(i_,j_)-1)-eps; -rho(i_,j_)*(rho(i_,j_)-1)-eps];
        g = [g; term(i_,j_)*(term(i_,j_)-1)-eps; -term(i_,j_)*(term(i_,j_)-1)-eps];
       
    end
end
g = [g; sum1(sum2(rho))-num_AMs; -sum1(sum2(rho))+num_AMs ];
g = [g; sum1(sum2(term)) - 1; -sum1(sum2(term)) + 1]; % One terminal


g = [g; X(:,1) - x_0; -X(:,1) + x_0];
g = [g; X(:,N+1) - x_f ; -X(:,N+1) + x_f];

%% Options
nlp_prob = struct('x', opt_variables, 'f', obj, 'g', g);
nlp_opts = struct;
nlp_opts.ipopt.print_level = 5;
nlp_opts.ipopt.tol = 1e-6;
nlp_opts.ipopt.max_iter =1000;
nlp_opts.ipopt.mu_strategy = 'adaptive';
nlp_opts.ipopt.linear_solver = 'mumps';
%nlp_opts.ipopt.jacobian_approximation = 'exact'; 
nlp_opts.ipopt.hessian_approximation = 'limited-memory';
nlp_opts.print_time = 1;

solver = nlpsol('solver', 'ipopt', nlp_prob, nlp_opts);
%% solve
LBx = [ ones(K*L*6,1)*-inf; ...
        repmat([ -0.1*pi; -pi/3; -pi/2 ; -pi/2; -inf; -inf; -inf; -inf], N+1, 1); ...
        repmat([tau_min; tau_min; u_min*ones(nu-2,1)], N, 1)];
UBx = [ ones(K*L*6,1)*inf; ...
        repmat([ pi/2; pi/3; pi/2 ; pi/2; inf; inf; inf; inf], N+1, 1); ...
        repmat([tau_max; tau_max; u_max*ones(nu-2,1)], N, 1)];

LBx = [ ones(K*L*6,1)*-inf; ...
        repmat([ -0.1*pi; -pi/3; -pi/2 ; -pi/2; -inf; -inf; -inf; -inf], N+1, 1); ...
        repmat([-inf; -inf; -inf*ones(nu-2,1)], N, 1)];
UBx = [ ones(K*L*6,1)*inf; ...
        repmat([ pi/2; pi/3; pi/2 ; pi/2; inf; inf; inf; inf], N+1, 1); ...
        repmat([inf; inf; inf*ones(nu-2,1)], N, 1)];

sol = solver('x0',  X_init_guess, ... 
             'lbx', -inf, 'ubx', +inf,...
             'lbg', -inf, 'ubg', 0);


% Extract the optimal solution
solution = full(sol.x);
rho_opt =reshape(solution(1:K*L), K, L);
flow_opt1 =reshape(solution(K*L+1:K*L*2), K, L);
flow_opt2 =reshape(solution(K*L*2+1:K*L*3), K, L);
flow_opt3 =reshape(solution(K*L*3+1:K*L*4), K, L);
flow_opt4 =reshape(solution(K*L*4+1:K*L*5), K, L);
term_opt = reshape(solution(K*L*5+1:K*L*6) ,K,L);
x_opt = reshape(solution(K*L*6+1 : K*L*6+nx*(N+1)), nx, N+1)';
u_opt = reshape(solution(( K*L*6+nx*(N+1)+1 ):end), nu, N)';
% Get the informations
optimal_value = full(sol.f);
exit_flag = solver.stats.success;
processing_time = solver.stats.t_wall_total;
disp(solver.stats.return_status)
end
%%
%Define Dynamic parameters and shapes
n = 4;
% dynamic parameters of each module 
params = define_params();
m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5}; thrust_limit= params{6};kt=params{7};c_1=params{8};c_2=params{9};

dh = [0,0,0.95,0;   % [alpha, a, d, theta]
      -pi/2, 0.9 , 0,0;
      0,-0.1,0.23,pi/2;
      pi/2,0,0,-pi/2;
      pi/2,0,0,0];
gravity = [0;0;-9.81];

% NLP parameters
dt = 0.1;
N = 100;

x_0 = [0;0;0;0;0;0;0;0];
x_f = [pi/4;0;0;0;0;0;0;0];
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

thrust_scale = 30;
u_max = thrust_limit *thrust_scale;
u_min = thrust_limit *(-thrust_scale);
tau_scale = 2.5;
tau_min = -0.2 *tau_scale ; 
tau_max =  0.2 *tau_scale ;

%initial guess
x_interp = zeros(N+1, 8);
vel_x1 = (x_f(1) - x_0(1))/(N*dt); vel_x2 = (x_f(2) - x_0(2))/(N*dt);
vel_x3 = (x_f(3) - x_0(3))/(N*dt); vel_x4 = (x_f(4) - x_0(4))/(N*dt);
for k = 1:8
    x_interp(:, k) = linspace(x_0(k), x_f(k), N+1)';
end
for k = 1:(N+1)
    x_interp(k, 5:8) = [vel_x1;vel_x2;vel_x3;vel_x4];
end

num_AMs = 4;
%K = 2*num_AMs+1; L = num_AMs; core = [num_AMs+1,1];
K=9;L=4;core=[5,1];
nu = 2 + K*L*4; zero_us = zeros(nu,1); 
rho_init = ones(K,L)/K/L*(num_AMs-1);
rho_init(core(1),core(2))=1;
X_init_guess = [reshape(rho_init,K*L,1);zeros(K*L*5,1);reshape(x_interp',(N+1)*8,1);repmat(zero_us, N, 1)];


%Solve NLP
[rho_opt, term_opt,flow_opt1,flow_opt2,flow_opt3,flow_opt4, x_opt, u_opt, optimal_value,exit_flag,processing_time] ....
          =  solve_NLP(params,num_AMs,K,L,core,dh,gravity,qo_desired,tau_min, tau_max, u_min,u_max,x_0,x_f,X_init_guess,dt,N);
disp(exit_flag)
disp(optimal_value)
disp(rho_opt)
%% plot
close all
figure('Position',[800,100,800,600])
time = (1:1:N+1)*dt;
subplot(2,1,1)
plot(time, x_opt(:,1:4))
hold on
plot(time, qo_desired, "--")
legend("q_1","q_2","q_3","q_4","q_{1,desired}","q_{2,desired}")
title("states")
axis tight;

subplot(2,1,2)
plot(time(1:end-1),u_opt(:,1:2))
legend("u_1","u_2")
title("motor inputs")
axis tight;
