function [lau_opt, x_opt, u_opt, optimal_value,exit_flag,processing_time, lau_history ,exit_flag_history, time_history] ....
          =  solve_NLP(params,num_AMs,K,L,core,dh,gravity,qo_desired,tau_min, tau_max, u_min,u_max,x_0,x_f,X_init_guess,dt,N ,max_iter,eps,gamma)
    addpath("../../casadi-3.6.7-windows64-matlab2018b")
    import casadi.*
    
    n = 4;
    m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5}; kt= params{7};c_1=params{8};c_2=params{9};
    mass_door = params{10}; handle_factor = params{11};

    nx = n*2;
    nu = 2 + K*L*4;
    
    x = MX.sym('x', nx, 1); % [q; qd];
    u = MX.sym('u', nu, 1);
    lau = MX.sym('lau', K,L); %shape
    
    % Dynamics
    [AM_com, AM_mass, AM_inertia] = get_inertia(lau,K,L, core ,m0, I0, d);
    mass =  {mass_door(1), mass_door(2), mass_door(3), AM_mass};
    inertia = {eye(3)*1, eye(3)*0.01, eye(3)*0.01, AM_inertia, zeros(3,3)};
    r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};
    
    tau = [-c_1*x(5);(-c_2*x(6) -kt*x(2) + mass{2} *handle_factor ); u(1); u(2)] ;
    F_ext = map_u2wrench( u(3:end), lau,K,L, core , mu , r , d);
    
    qdd = FD(n, dh, mass, inertia, r_i_ci, gravity, x(1:4), x(5:8), tau, F_ext);
    
    x_dot = [x(5:8);qdd];
    x_next = x + dt*x_dot ; 
    
    % Create CasADi function for the dynamics
    f = Function('f', {x, u, lau}, {x_next});
    
    U = MX.sym('U', nu, N);
    X = MX.sym('X', nx, N+1);
    opt_variables = [reshape(lau, K*L,1);reshape(X, nx*(N+1), 1); reshape(U, nu*N, 1)];
    %%
    success=0;
    iter = 1;
    
    lau_history = cell(1,100);
    exit_flag_history = cell(1,100);
    time_history = cell(1,100);
    obj_history = cell(1,100);
    
    tic
    while   ( success==0 || eps> (0.005*gamma) ) && iter < 15 
        fprintf("num AM: %d, iter: %d, eps: %f\n",num_AMs, iter,eps);
        obj = 0;
        g = [];
        Q1 = diag([0,0,0,0,1,1,1,1])*0.1; Q2 =diag([1,1]); Q3 = diag([1,1,1,1,1,1,1,1])*1e1;
        R = diag(ones(nu,1))*1;
        
        for k = 1:N
            if k <= 1/dt
              %g = [g; X(:,k) - x_0 ; -X(:,k) + x_0]; 
              obj = obj + (X(:,k) - x_0)'*Q3*(X(:,k) - x_0);
            elseif N-k <= 1/dt
              %g = [g; X(:,k) - x_f ; -X(:,k) + x_f];
              obj = obj + (X(:,k) - x_f)'*Q3*(X(:,k) - x_f);
            else
              obj = obj + X(:,k)'*Q1*X(:,k) + (X(1:2,k) - qo_desired(:,k))'*Q2*(X(1:2,k) - qo_desired(:,k));
              obj = obj + U(:,k)'*R*U(:,k);
            end
            % dynamics equality constraint
            g = [g; X(:,k+1) - f(X(:,k), U(:,k),lau) ; -X(:,k+1) + f(X(:,k), U(:,k),lau)]; 
        end
        
        for i_=1:K
            for j_=1:L
                if(i_==core(1) && j_==core(2))
                    g = [g; lau(i_,j_)-1; -lau(i_,j_)+1];
                else
                    %neighbor_sum = MX(0); 
                    neighbor_max = MX(0);
                    if i_ > 1 && i_>core(1)
                        %neighbor_sum = neighbor_sum + lau(i_-1, j_);
                        neighbor_max = fmax(neighbor_max, lau(i_-1, j_));
                    end
                    if i_ < K && i_<core(1)
                        %neighbor_sum = neighbor_sum + lau(i_+1, j_);
                        neighbor_max = fmax(neighbor_max, lau(i_+1, j_));
                    end
                    if j_ > 1 && j_>core(2)
                        %neighbor_sum = neighbor_sum + lau(i_, j_-1);
                        neighbor_max = fmax(neighbor_max, lau(i_, j_-1));
                    end
                    if j_ < L && j_<core(2)
                        %neighbor_sum = neighbor_sum + lau(i_, j_+1);
                        neighbor_max = fmax(neighbor_max, lau(i_, j_+1));
                    end
                    
                    %g = [g; -lau(i_,j_); lau(i_,j_) - fmin(MX(1),neighbor_sum + eps) ];
                    g = [g; -lau(i_,j_); lau(i_,j_) - neighbor_max ];
                end
                
                g = [g; lau(i_,j_)*(lau(i_,j_)-1)-eps; -lau(i_,j_)*(lau(i_,j_)-1)-eps];
            end
        end
        g = [g; sum1(sum2(lau))-num_AMs; -sum1(sum2(lau))+num_AMs ];
        %g = [g; sum1(sum2(lau))-num_AMs];
        
        g = [g; X(:,1) - x_0; -X(:,1) + x_0];
        g = [g; X(:,N+1) - x_f ; -X(:,N+1) + x_f];
        
        % solver options
        nlp_prob = struct('x', opt_variables, 'f', obj, 'g', g);
        nlp_opts = struct;
        nlp_opts.ipopt.print_level = 1;
        nlp_opts.ipopt.tol = 1e-3;
        nlp_opts.ipopt.max_iter = max_iter;%*(1+0.15*iter);
        nlp_opts.ipopt.mu_strategy = 'monotone';
        nlp_opts.ipopt.linear_solver = 'mumps';
        %nlp_opts.ipopt.jacobian_approximation = 'exact'; 
        nlp_opts.ipopt.hessian_approximation = 'limited-memory';
        nlp_opts.print_time = 1;
        
        solver = nlpsol('solver', 'ipopt', nlp_prob, nlp_opts);
        
        % solve
        LBx = [ ones(K*L,1)*-inf; ...
                repmat([ -0.1*pi; -pi/3; -pi/2 ; -pi/2; -inf; -inf; -inf; -inf], N+1, 1); ...
                repmat([tau_min; tau_min; u_min*ones(nu-2,1)], N, 1)];
        UBx = [ ones(K*L,1)*inf; ...
                repmat([ pi/2; pi/3; pi/2 ; pi/2; inf; inf; inf; inf], N+1, 1); ...
                repmat([tau_max; tau_max; u_max*ones(nu-2,1)], N, 1)];
        
        if iter==1
            sol = solver('x0',  X_init_guess, ... 
                         'lbx', LBx, 'ubx', UBx,...
                         'lbg', -inf, 'ubg', 0);
        else
            sol = solver('x0',  sol.x, ... 
                         'lbx', LBx, 'ubx', UBx,...
                         'lbg', -inf, 'ubg', 0);
        end
        
        disp(solver.stats.return_status)
        full_sol = full(sol.x);
        current_lau = reshape(full_sol(1:K*L), K, L);
        lau_history{iter} = current_lau;
        exit_flag_history{iter} = solver.stats.return_status;
        time_history{iter} = solver.stats.t_wall_total;
        obj_history{iter} = full(sol.f);
        fprintf("current object function value: %f\n", full(sol.f));
        disp(current_lau)
    
        success = solver.stats.success;
        eps = eps*gamma;
        iter= iter+1;
    end
    toc
    
    solution = full(sol.x);
    lau_opt =reshape(solution(1:K*L), K, L);
    x_opt = reshape(solution(K*L+1 : K*L+nx*(N+1)), nx, N+1)';
    u_opt = reshape(solution(( K*L+nx*(N+1)+1 ):end), nu, N)';
    optimal_value = full(sol.f);
    exit_flag = solver.stats.success;
    processing_time = sum(cell2mat(time_history));
    disp(solver.stats.return_status)
end

%%
% Define Dynamic parameters and shapes
n = 4;
params = define_params();
m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5}; thrust_limit= params{6};kt=params{7};c_1=params{8};c_2=params{9};
mass_door = params{10}; handle_factor = params{11};
dh = [0,0,0.95,0;   % [alpha, a, d, theta]
      -pi/2, 0.9 , 0,0;
      0,-0.1,0.23,pi/2;
      pi/2,0,0,-pi/2;
      pi/2,0,0,0];
gravity = [0;0;-9.81];

% NLP parameters
dt = 0.1;
N = 80;

x_0 = [0;0;0;0;0;0;0;0];
x_f = [pi/4;0;0;0;0;0;0;0];
qo_desired = zeros(2,N+1);
T = N*dt; t0 = 1; t1 = 1; %1sec
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
qo_desired = [repmat(x_0(1:2),1,t1/dt), qo_desired, repmat(x_f(1:2),1,t1/dt)];
N = N + t1/dt*2;

thrust_scale = 3;
tau_scale = 0;

u_max = thrust_limit *thrust_scale;
u_min = thrust_limit *(-thrust_scale);
tau_min = -0.2 *tau_scale ; 
tau_max =  0.2 *tau_scale ;

% initial guess
x_interp = zeros(N+1, 8);
vel_x1 = (x_f(1) - x_0(1))/(N*dt); vel_x2 = (x_f(2) - x_0(2))/(N*dt);
vel_x3 = (x_f(3) - x_0(3))/(N*dt); vel_x4 = (x_f(4) - x_0(4))/(N*dt);
for k = 1:8
    x_interp(:, k) = linspace(x_0(k), x_f(k), N+1)';
end
for k = 1:(N+1)
    x_interp(k, 5:8) = [vel_x1;vel_x2;vel_x3;vel_x4];
end
for MAX_ITER=[1000]
    fprintf("\nmax_iter: %d\n", MAX_ITER);
        for num_AMs = 11:1:15
            %num_AMs = 12;
            K = 2*num_AMs-1; L = num_AMs; core = [num_AMs,1];
            %K=9; L=5; core=[5,1];
            nu = 2 + K*L*4; zero_us = zeros(nu,1); 
            lau_init = ones(K,L)/K/L*(num_AMs-1);
            lau_init(core(1),core(2))=1; %lau_init(core(1)+1,core(2))=1.0; %lau_init(core(1)+1,core(2)+1)=0.5; 
            X_init_guess = [reshape(lau_init,K*L,1);reshape(x_interp',(N+1)*8,1);repmat(zero_us, N, 1)];
            
            max_iter = MAX_ITER;%default 50
            eps = 0.25;
            gamma = 0.3;
        
            %Solve NLP
            [lau_opt, x_opt, u_opt, optimal_value, exit_flag, processing_time, lau_history, exit_flag_history, time_history] ....
                      =  solve_NLP(params,num_AMs,K,L,core,dh,gravity,qo_desired,tau_min, tau_max, u_min,u_max,x_0,x_f,X_init_guess,dt,N,max_iter,eps,gamma);
            fprintf("num AMs : %d\n", num_AMs);
            fprintf("exit flag: %d \n", exit_flag);
            fprintf("optimal value: %f \n", optimal_value);
            fprintf("lau: \n"); disp(lau_opt)
            
            filename = sprintf('result/hovor/max_iter_%d/%d_%d_%d.mat', max_iter, num_AMs, thrust_scale , tau_scale);
            save(filename);
        end
end
%% plot
close all
figure('Position',[900,100,900,800])
time = (1:1:N+1)*dt;
subplot(4,1,1)
plot(time, x_opt(:,1:4))
hold on
plot(time, qo_desired, "--")
legend("q_1","q_2","q_3","q_4","q_{1,ref}","q_{2,ref}")
title("states")
axis tight;

subplot(4,1,2)
plot(time(1:end-1),u_opt(:,1:2))
legend("u_1","u_2")
title("motor inputs")
axis tight;

[AM_com, AM_mass, AM_inertia]  = get_inertia_double(lau_opt,K,L, core ,m0, I0, d);
mass =  {mass_door(1), mass_door(2), mass_door(1), AM_mass};
inertia = {eye(3)*1, eye(3)*0.1, eye(3)*0.1, AM_inertia, zeros(3,3)};
r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};

wrench = zeros(N,6);
tau = zeros(N,n);
for i=1:N
    wrench(i,:) = map_u2wrench_double( u_opt(i,3:end)',lau_opt,K,L, core , mu , r , d);
    q = x_opt(i,1:4)'; qd = x_opt(i,5:8)'; qdd = (x_opt(i+1,5:8) -x_opt(i+1,5:8) )'/dt; 
    F_ext = wrench(i,:)';
    tau(i,:) =  [-c_1*qd(1);(-c_2*qd(2) -kt*q(2) + mass{2}*handle_factor); u_opt(i,1); u_opt(i,2)] + ...
        newton_euler_inverse_dynamics_double(n, dh, mass, inertia, r_i_ci, gravity, q, qd, qdd, F_ext);
end
subplot(4,1,3)
plot(time(1:end-1),wrench)
legend("m_x","m_y","m_z","f_x","f_y","f_z")
title("Wrench exp. AM frame")
axis tight;

subplot(4,1,4)
plot(time(1:end-1),tau)
axis tight
legend("tau1","tau2","tau3","tau4")
title("generalized force")
axis tight;

%plot 3d video
do_view=1; q =  [0;0;0;0]; g=[0;0;-9.81];
robot = generate_door(n,dh,r_i_ci,d, g, lau_opt, core, mass,inertia, do_view,q);

slow_factor =1; force_scale = 0.2;
%save_plot_tree(robot,dh, params, x_opt,u_opt, dt,N,slow_factor, force_scale, lau_opt, core, K, L)
plot_tree(robot,dh, params, x_opt,u_opt, dt,N,slow_factor, force_scale, lau_opt, core, K, L)
