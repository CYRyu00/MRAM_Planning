addpath("../../casadi-3.6.7-windows64-matlab2018b" , "dynamics", "params")
%%
load("..\3D_ver2\data\result_9_5\hover\max_iter_1000\5_3_0.mat")
shape = zeros([K,L]);
shape_idx = zeros([K,L]);
x_d = x_opt;
u_d = zeros([N,4*num_AMs]);
[AMs_rows, AMs_cols] = find(lau_opt >= 0.9);
for i = 1:length(AMs_rows)
    shape(AMs_rows(i),AMs_cols(i)) = 1;
    shape_idx(AMs_rows(i),AMs_cols(i)) = i;
    u_d(:,4*i-3:4*i) = u_opt(:, 4*((AMs_rows(i)-1)*L + AMs_cols(i))-1 : 4*((AMs_rows(i)-1)*L + AMs_cols(i))+2 );
end
shape(core(1),core(2)) = 2;
fprintf('Shape : \n')
disp(shape)

%% CASDI function: Get A and B
import casadi.*

delta_inertia = MX.sym('delta_inertia',1,1);
delta_k = MX.sym('delta_k',1,1);
disturb = MX.sym('disturb',n,1);
params = define_params();
m0 = params{1} *delta_inertia; I0 = params{2}*delta_inertia; mass_door = params{10} *delta_inertia;
mu = params{3}; r= params{4}; d= params{5}; thrust_limit= params{6};kt=params{7} *delta_k;c_1=params{8};c_2=params{9};

[AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);
mass =  {mass_door(1), mass_door(2), mass_door(3), AM_mass};
inertia = {eye(3)*1, eye(3)*0.1, eye(3)*0.1, AM_inertia, zeros(3,3)};
r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};

nx = n*2;
nu = num_AMs*4;
x = MX.sym('x', nx, 1); % [q; qd];
u = MX.sym('u', nu, 1);
tau = [-c_1*x(5);(-c_2*x(6) -kt*x(2) + mass{2}*handle_factor); 0; 0] + disturb;
F_ext = map_u2wrench( u, shape , mu , r , d);

qdd = FD_ver2(n, dh, mass, inertia, r_i_ci, gravity, x(1:4), x(5:8), tau, F_ext);

x_dot = [x(5:8);qdd];

A = jacobian(x_dot, x);
B = jacobian(x_dot, u);

x_dot_func = Function('x_dot_func', {x ,u ,delta_inertia, delta_k, disturb }, {x_dot});
A_func = Function('A_func', {x, u, delta_inertia ,delta_k,disturb }, {A});
B_func = Function('B_func', {x, u, delta_inertia, delta_k,disturb }, {B});
%%
%hover 1sec at first and last
%params = define_params();
%m0 = params{1}; I0 = params{2}; mass_door = params{10};
%[AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);
%N = N + 1/dt*2;
%x_d = [repmat(x_0', 1/dt,1); x_d ; repmat(x_f', 1/dt,1)];
%u_hovor = ones(1,nu)*(AM_mass + mass_door(3))*norm(gravity)/num_AMs/4; 
%u_d = [repmat(u_hovor,1/dt,1 ); u_d ; repmat(u_hovor,1/dt,1 )];

dt_sim = 0.01;
N_sim = N*dt/dt_sim;

t = linspace(0, dt*N, N+1);         
t_sim = linspace(0, dt_sim*N_sim, N_sim+1);  

% Interpolation
x_d_interp = interp1(t, x_d, t_sim, 'pchip');  % 1001x8
u_d_interp = interp1(t(1:end-1), u_d, t_sim(1:end-1), 'pchip');  % 1000x40

close all
figure('Position',[100 500 600 400])
subplot(2,1,1)
plot(t, x_d,'.');
title("x_d, dt = 0.1")
subplot(2,1,2)
plot(t_sim, x_d_interp,'.','MarkerSize',2);
title("dt = 0.01")

figure('Position',[700 500 600 400])
subplot(2,1,1)
plot(t(1:end-1), u_d(:,1:8),'.');
title("u_d, dt = 0.1")
subplot(2,1,2)
plot(t_sim(1:end-1), u_d_interp(:,1:8),'.','MarkerSize',2);
title("dt = 0.01")

%%
fprintf("Check Controllability for all x_d, u_d \n")
A_arr = cell(1,N_sim);
B_arr = cell(1,N_sim);
delta_inertia = 1.0;
delta_k = 1.0;
disturb = [0;0;0;0];
for i=1:N_sim
    A_val = full(A_func(x_d_interp(i,:), u_d_interp(i,:),delta_inertia, delta_k,disturb ));  % nx x nx
    B_val = full(B_func(x_d_interp(i,:), u_d_interp(i,:),delta_inertia, delta_k,disturb ));  % nx x nu
    A_arr{i} = A_val; B_arr{i} = B_val;
    C = ctrb(A_val, B_val);
    rank_C = rank(C);
    is_controllable = (rank_C == size(A_val,1));
    if is_controllable == false
        break
    end
end
if is_controllable
    disp('✅ System is controllable.');
else
    fprintf('❌ System is NOT controllable. at time step : %d\n\n',i)
end
%%
K_arr = cell(1, N_sim);
Q = diag(ones(nx,1));
Qf = diag(ones(nx,1));
R = diag(ones(nu,1))*1e-2;
P = Qf;  % Terminal cost

Ad_arr = cell(1, N_sim);
Bd_arr = cell(1, N_sim);
for i = 1:N_sim
    A = A_arr{i}; B = B_arr{i};
    M = [A, B; zeros(size(B,2), size(A,1)+size(B,2))];
    Md = expm(M * dt);  % dt: discretization time step
    Ad_arr{i} = Md(1:size(A,1), 1:size(A,2));
    Bd_arr{i} = Md(1:size(A,1), size(A,2)+1:end);
end

for i = N_sim:-1:1
    A = Ad_arr{i};
    B = Bd_arr{i};
    S = R + B' * P * B;
    K = (S) \ (B' * P * A);  % Optimal gain
    K_arr{i} = K;
    P = Q + A' * P * A - A' * P * B * (S \ (B' * P * A));
end
%% Simulation 
%TODO
x_sim = zeros(N_sim+1,nx);
u_sim = zeros(N_sim,nu);
x_sim(1,:)=x_d_interp(1,:);
disturb_sim = zeros(N_sim, n);
delta_inertia = 1.0;
delta_k = 1.0;
sigma = 0.10;
mean = 0.0;
max_val = 0.05;
disturb = mean*ones(n,1);

rng('shuffle') 
for i = 1:N_sim
    %if mod(i,10) == 1
        disturb_dot = randn(n, 1) *sigma*3;
        disturb = disturb + disturb_dot *dt_sim;
        
        disturb = min(max(disturb, -max_val), max_val) ;
    %end
    
    disturb_sim(i,:) = disturb;
    u_sim(i,:) = ( u_d_interp(i,:)' - K_arr{i}*(x_sim(i,:)' - x_d_interp(i,:)') )';
    x_sim(i+1,:) = x_sim(i,:) + full(x_dot_func(x_sim(i,:),u_sim(i,:),delta_inertia, delta_k, disturb))'*dt_sim;
end

% Plot
%figure('Position',[100 100 1000 800])
figure('Position',[100 100 800 600])
colors = lines(nx);
legend_entries = cell(1, nx/2);
subplot(3,1,1)
hold on
for i = 1:nx/2
    % Plot simulation
    plot(t_sim, x_sim(:,i), 'Color', colors(i,:), 'LineWidth', 1.0);
    legend_entries{2*i-1} = sprintf('$q_%d$', i);
    
    % Plot reference
    plot(t_sim, x_d_interp(:,i), '--', 'Color', colors(i,:), 'LineWidth', 1.0);
    legend_entries{2*i} = sprintf('$q_{%d,d}$', i);
end
legend(legend_entries, 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 10)
xlabel('Time [s]'); ylabel('State Value')
title_name = sprintf(['x_{sim} vs x_{desired}, \nMass X %.2f, Spring x %.2f,\n Disturbance : ' ...
                        'Sigma: %.2f, Mean: %.2f, Max val :%.2f'],delta_inertia, delta_k , sigma, mean ,max_val);
title(title_name , 'FontSize',10)
grid on

subplot(3,1,2)
legend_entries = cell(1, nx/2);
hold on
for i = (nx/2+1):nx
    % Plot simulation
    plot(t_sim, x_sim(:,i), 'Color', colors(i,:), 'LineWidth', 1.0);
    legend_entries{2*i-1 -nx} = sprintf('$\\dot{q}_%d$', i);
    
    % Plot reference
    plot(t_sim, x_d_interp(:,i), '--', 'Color', colors(i,:), 'LineWidth', 1.0);
    legend_entries{2*i -nx} = sprintf('$\\dot{q}_{%d,d}$', i);
end
legend(legend_entries, 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 10)
xlabel('Time [s]'); ylabel('State Value')
grid on

subplot(3,1,3)
colors = lines(nu);
hold on
for i = 1:nu
    plot(t_sim(1:end-1), u_sim(:,i), 'Color', colors(i,:), 'LineWidth', 1.0);    
    plot(t_sim(1:end-1), u_d_interp(:,i), '--', 'Color', colors(i,:), 'LineWidth', 1.0);
end
title("u_{sim} vs u_{desired}")
xlabel('Time [s]'); ylabel('Input Value')
grid on

wrench = zeros(N_sim,6);
tau = zeros(N_sim,n);
params = define_params();
m0 = params{1}; I0 = params{2}; kt_double = params{7} *delta_k; mass_door = params{10} *delta_inertia; 
[AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);

mass_double =  {mass_door(1), mass_door(2), mass_door(3), AM_mass};
inertia_double = {eye(3)*1, eye(3)*0.1, eye(3)*0.1, AM_inertia, zeros(3,3)};
r_i_ci_double = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};

for i=1:N_sim
    wrench(i,:) = map_u2wrench_double( u_sim(i,:)', shape , mu , r , d);
    q = x_sim(i,1:4)'; qd = x_sim(i,5:8)'; qdd = (x_sim(i+1,5:8) -x_sim(i+1,5:8) )'/dt; 
    F_ext = wrench(i,:)';
    tau(i,:) =  [-c_1*qd(1); full(-c_2*qd(2) -kt_double*q(2) + mass_double{2}*handle_factor); u_sim(i,1); u_sim(i,2)] + ...
        newton_euler_inverse_dynamics_double(n, dh, mass_double, inertia_double, r_i_ci_double, gravity, q, qd, qdd, F_ext) ...
            + disturb_sim(i,:)';
end
%figure('Position',[1100 100 700 600])
figure('Position',[1100 100 800 600])
subplot(2,1,1)
plot(t_sim(1:end-1),wrench)

legend({'$m_x$', '$m_y$', '$m_z$', '$f_x$', '$f_y$', '$f_z$'}, ...
       'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 10);
title("Wrench exp. AM frame")
axis tight;
grid on

subplot(2,1,2)
colors = lines(n);
hold on
for i = 1:n
    plot(t_sim(1:end-1), tau(:,i), 'Color', colors(i,:), 'LineWidth', 1.0);    
    plot(t_sim(1:end-1), disturb_sim(:,i), '.', 'Color', colors(i,:), 'LineWidth', 1.0);
end
axis tight
legend({'$\tau_1$','$d_1$', '$\tau_2$','$d_2$', '$\tau_3$', '$d_3$','$\tau_4$','$d_4$',}, ...
       'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 10);

title("generalized force")
axis tight;
grid on
