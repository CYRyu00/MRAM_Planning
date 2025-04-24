addpath("dynamics","../params" ,"plot")
clear
% Define Dynamic parameters and shapes
n = 4;
params = define_params();
m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5}; thrust_limit= params{6};kt=params{7};c_1=params{8};c_2=params{9}; mass_door = params{10}; handle_factor = params{11};
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
qo_desired = zeros(4,N+1);
T = N*dt; t0 = 1; t1 = 1; %1sec
for j=1:N+1
    t = (j-1)*dt;
    if t < t0 
        qo_desired(:,j) = [0 ; pi/6*(cos(pi/t0*t) - 1)/2; 0; pi/4*(cos(pi/t0*t) - 1)/2 ];
    elseif t < T-t0
        qo_desired(:,j) = [(x_f(1)-x_0(1))/(T-t0)*(t-t0) ;-pi/6;0 ; -pi/4];
    else
        qo_desired(:,j) = [(x_f(1)-x_0(1))/(T-t0)*(t-t0) ;pi/6*(cos(pi/t0*t -pi/t0*T) - 1)/2; ...
                            0; pi/4*(cos(pi/t0*t -pi/t0*T) - 1)/2];
    end
end
qo_desired = [repmat(x_0(1:4),1,t1/dt), qo_desired, repmat(x_f(1:4),1,t1/dt)];
N = N + t1/dt*2;

% initial guess
x_interp = zeros(N+1, 8);
vel_x1 = (x_f(1) - x_0(1))/( (N-2*t1)*dt); vel_x2 = (x_f(2) - x_0(2))/( (N-2*t1)*dt);
vel_x3 = (x_f(3) - x_0(3))/( (N-2*t1)*dt); vel_x4 = (x_f(4) - x_0(4))/( (N-2*t1)*dt);
for k = 1:8
    x_interp(t1/dt:(end-1-t1/dt), k) = linspace(x_0(k), x_f(k), N+1-2*t1/dt)';
end
for k = t1/dt:(N-t1/dt)
    x_interp(k, 5:8) = [vel_x1;vel_x2;vel_x3;vel_x4];
end

CASE = 1;
thrust_scale = 1;
tau_scale = 0;

u_max = thrust_limit *thrust_scale;
u_min = thrust_limit *(-thrust_scale);
tau_min = -0.2 *tau_scale ; 
tau_max =  0.2 *tau_scale ;

max_iter = 2000;
eps = 0.25;
gamma = 0.3;

fprintf("\nmax_iter: %d\n", max_iter);
    for num_AMs = 3:1:3
        if CASE == 1
            K = 2*num_AMs-1; L = num_AMs; core = [num_AMs,1];
            filename = sprintf('data/test.mat');
            %mkdir data/result/hover/new_ref
            %filename = sprintf('data/result/hover/new_ref/%d_%d_%d.mat', num_AMs, thrust_scale , tau_scale);
            %filename = sprintf('data/result/hover/max_iter_%d/%d_%d_%d.mat', max_iter, num_AMs, thrust_scale , tau_scale);
        elseif CASE == 2
            K=9; L=5 ; core=[5,1];
            mkdir data/result_9_5/hover/new_ref_2
            filename = sprintf('data/result_9_5/hover/new_ref_2/%d_%d_%d.mat', num_AMs, thrust_scale , tau_scale);
            %filename = sprintf('data/result_9_5/hover/max_iter_%d/%d_%d_%d.mat', max_iter, num_AMs, thrust_scale , tau_scale);
        elseif CASE == 3
            K=11; L=6 ; core=[6,1];
            filename = sprintf('data/result_11_6/hover/max_iter_%d/%d_%d_%d.mat', max_iter, num_AMs, thrust_scale , tau_scale);
        end
        fprintf("[K, L] = [%d , %d]\n", K,L);
        
        nu = 2 + K*L*4; zero_us = zeros(nu,1); 
        rho_init = ones(K,L)/K/L*(num_AMs-1);
        rho_init(core(1),core(2))=1;
        X_init_guess = [reshape(rho_init,K*L,1);reshape(x_interp',(N+1)*8,1);repmat(zero_us, N, 1)];
        
        %Solve NLP
        [rho_opt, x_opt, u_opt, optimal_value, exit_flag, processing_time, rho_history, exit_flag_history, time_history] ....
                  =  solve_nlp(params,num_AMs,K,L,core,dh,gravity,qo_desired,tau_min, tau_max, u_min,u_max,x_0,x_f,X_init_guess,dt,N,max_iter,eps,gamma);
        fprintf("num AMs : %d\n", num_AMs);
        fprintf("exit flag: %d \n", exit_flag);
        fprintf("optimal value: %f \n", optimal_value);
        fprintf("rho: \n"); disp(rho_opt)
       
        save(filename);
    end
%% plot
plot_graph_video
