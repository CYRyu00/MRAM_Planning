addpath("../../casadi-3.6.7-windows64-matlab2018b" , "dynamics", "params")
num_AMs_arr = 1:5;
min_eigval_arr = []; 
xT_W_r_x_cell = cell(1, length(num_AMs_arr));

for num_AMs = num_AMs_arr
    filename = sprintf('../3D_ver2/data/result_9_5/hover/c1_10_times/%d_1_0.mat', num_AMs);
    load(filename)
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
    fprintf('num AMs: %d, Shape : \n', num_AMs)
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
    dt_sim = 0.01;
    N_sim = N*dt/dt_sim;
    
    t = linspace(0, dt*N, N+1);         
    t_sim = linspace(0, dt_sim*N_sim, N_sim+1);  
    
    x_d_interp = interp1(t, x_d, t_sim, 'pchip');  
    u_d_interp = interp1(t(1:end-1), u_d, t_sim(1:end-1), 'pchip');  
    
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
    N_horizon = 10;  % Set reachability window

    for i = 1:N_sim - N_horizon + 1
        W_r = zeros(2*n);
        Phi = eye(2*n);
 
        for j = 0:N_horizon-1
            A_j = A_arr{i + j};
            B_j = B_arr{i + j};
    
            if j > 0
                Phi = A_arr{i + j - 1} * Phi;
            end
    
            W_r = W_r + Phi * B_j * B_j' * Phi';
        end
   
        [V, D] = eig(W_r); 
        [min_eigval, min_idx] = min(diag(D)); 
        min_eigval_arr(i, num_AMs) = min_eigval;

        xT_W_r_x = zeros(1,2*n);
        for k = 1:2*n
            e_k = zeros(2*n,1); e_k(k) = 1;
            xT_W_r_x(k) = e_k'*W_r*e_k;
        end
        xT_W_r_x_cell{num_AMs} = [xT_W_r_x_cell{num_AMs}; xT_W_r_x];

    end
end
%% plot
close all
t_tmp = (1:N_sim - N_horizon + 1)*dt_sim;
legend_entries = arrayfun(@(k) ['$' num2str(k) '$'], num_AMs_arr, 'UniformOutput', false);

figure
plot(t_tmp, log(min_eigval_arr)/log(10))
legend(legend_entries, 'FontSize', 10, 'Interpreter', 'latex');
title('minimum eigenvalue for each number of AM' ...
    , 'FontSize', 14, 'Interpreter', 'latex');
xlabel('time[sec]', 'FontSize', 14, 'Interpreter', 'latex')
ylabel( '$\log_{10}{(\min{\lambda})}$','FontSize', 14,'Interpreter', 'latex')
grid on

figure 
hold on
for num_AMs = num_AMs_arr
    plot(t_tmp, log(xT_W_r_x_cell{num_AMs}(:,1))/log(10))
end
legend(legend_entries, 'FontSize', 10, 'Interpreter', 'latex');
title({'$e_1^T W_r e_1$ for each number of AM'}, ...
       'Interpreter','latex','FontSize',14);
xlabel('time[sec]', 'FontSize', 14, 'Interpreter', 'latex')
ylabel( '$\log_{10}{(e_1^T W_r e_1)}$','FontSize', 14,'Interpreter', 'latex')
grid on