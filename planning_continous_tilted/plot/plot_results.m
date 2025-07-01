close all
num_AMs_arr=[];
opt_val_arr=[];
time_arr =[];
for num_AMs = 1:1:3
    filename = sprintf('data/Q2_1e0_1110/10sec/%d_1_0.mat', num_AMs);
    load(filename);
    fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, processing_time)
   
    num_AMs_arr = [num_AMs_arr ; num_AMs];
    opt_val_arr = [opt_val_arr; optimal_value];
    time_arr = [time_arr; processing_time];
    
    %disp(rho_opt)
    do_view=1; g=[0;0;-9.81];q=[0;0;0;0];
    [AM_com, AM_mass, AM_inertia]  = get_inertia_double(rho_opt,K,L, core ,m0, I0, d);
    mass =  {10, 1, 0.5, AM_mass};
    inertia = {eye(3)*1, eye(3)*0.1, eye(3)*0.1, AM_inertia, zeros(3,3)};
    r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};
    robot = generate_door(n,dh,r_i_ci,d, g, rho_opt, core, mass,inertia, do_view,q);
end

figure
hold on
plot(num_AMs_arr, opt_val_arr,'o-');
plot(num_AMs_arr, time_arr/60,'o-');
grid on
hold off
xlabel("Number of Modules");
legend("f*", "time[min]")
%% old version
close all
for num_AMs = 1:1:5
    filename = sprintf('data/old_result/hover/max_iter_2000/%d_3_0.mat', num_AMs);
    load(filename);
    [optimal_value, argmin]  = min(cell2mat(all_optimal_value));
    fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, elapsed_time)
    shape = shapes{argmin};

    params = define_params();
    m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5}; thrust_limit= params{6};kt=params{7};c_1=params{8};c_2=params{9};
    
    [AM_com, AM_mass, AM_inertia] = get_inertia_old(shape ,m0, I0, d);
    mass =  {10, 1, 0.5, AM_mass};
    inertia = {eye(3)*1, eye(3)*0.1, eye(3)*0.1, AM_inertia, zeros(3,3)};
    r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};

    do_view=1; g=[0;0;-9.81];q=[0;0;0;0];n=4;
    dh = [0,0,0.95,0;   % [alpha, a, d, theta]
      -pi/2, 0.9 , 0,0;
      0,-0.1,0.23,pi/2;
      pi/2,0,0,-pi/2;
      pi/2,0,0,0];
    robot = generate_door_old(n,dh,r_i_ci, d, g, shape, mass,inertia, do_view,q);
end
%%
% optimal value
close all
plot_old_arr = 1:1:7;
plot_full_arr = 1:1:10;
plot_9_5_arr = 1:1:10;
plot_11_6_arr = 1:1:10;

num_AMs_arr=[];
opt_val_arr=[];
figure('Position',[300 500 500 400])

for num_AMs = plot_old_arr
    filename = sprintf('data/old_result/ref_1/%d_1_0.mat', num_AMs);
    load(filename);
    for idx = 1:length(cell2mat(all_optimal_value))
        if all_exit_flag{idx} == 0
            all_optimal_value{idx} = inf;
        end
    end
    optimal_value = min(cell2mat(all_optimal_value));
    fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, elapsed_time)
    
    num_AMs_arr = [num_AMs_arr ; num_AMs];
    opt_val_arr = [opt_val_arr; optimal_value];
end
hold on
plot(num_AMs_arr, opt_val_arr,'o-');

num_AMs_arr=[];
opt_val_arr=[];

for num_AMs = plot_full_arr
    filename = sprintf('data/result_full/ref_1/%d_1_0.mat', num_AMs);
    load(filename);
    fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, processing_time)
    
    num_AMs_arr = [num_AMs_arr ; num_AMs];
    opt_val_arr = [opt_val_arr; optimal_value];
end
plot(num_AMs_arr, opt_val_arr,'o-');

num_AMs_arr=[];
opt_val_arr=[];
for num_AMs = plot_9_5_arr
    filename = sprintf('data/result_9_5/ref_1/%d_1_0.mat', num_AMs);
    load(filename);
    fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, processing_time)
    num_AMs_arr = [num_AMs_arr ; num_AMs];
    opt_val_arr = [opt_val_arr; optimal_value];
end
plot(num_AMs_arr, opt_val_arr,'o-');

num_AMs_arr=[];
opt_val_arr=[];
for num_AMs = plot_11_6_arr
    filename = sprintf('data/result_11_6/ref_1/%d_1_0.mat', num_AMs);
    load(filename);
    fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, processing_time)
    
    num_AMs_arr = [num_AMs_arr ; num_AMs];
    opt_val_arr = [opt_val_arr; optimal_value];
end
plot(num_AMs_arr, opt_val_arr,'o-');

grid on
axis tight
hold off
xlabel("Number of Modules")
title("Object Function value : f*")
legend("Prev. Full", "New Full","New (9,5)","New (11,6)",'FontSize', 12 )

% ============================== Time ===================================
num_AMs_arr=[];
time_arr =[];
figure('Position',[800 500 500 400])

for num_AMs = plot_old_arr
    filename = sprintf('data/old_result/ref_1/%d_1_0.mat', num_AMs);
    load(filename);
    fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, elapsed_time)
    
    num_AMs_arr = [num_AMs_arr ; num_AMs];
    time_arr = [time_arr; elapsed_time];
end
hold on
plot(num_AMs_arr, time_arr/60,'o-');

num_AMs_arr=[];
time_arr =[];
for num_AMs = plot_full_arr
    filename = sprintf('data/result_full/ref_1/%d_1_0.mat', num_AMs);
    load(filename);
    fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, processing_time)
    
    num_AMs_arr = [num_AMs_arr ; num_AMs];
    time_arr = [time_arr; processing_time];
end
plot(num_AMs_arr,time_arr/60,'o-');

num_AMs_arr=[];
time_arr =[];
for num_AMs = plot_9_5_arr
    filename = sprintf('data/result_9_5/ref_1/%d_1_0.mat', num_AMs);
    load(filename);
    fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, processing_time)
    
    num_AMs_arr = [num_AMs_arr ; num_AMs];
    time_arr = [time_arr; processing_time];
end
plot(num_AMs_arr, time_arr/60,'o-');

num_AMs_arr=[];
time_arr =[];
for num_AMs = plot_11_6_arr
    filename = sprintf('data/result_11_6/ref_1/%d_1_0.mat', num_AMs);
    load(filename);
    fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, processing_time)
    
    num_AMs_arr = [num_AMs_arr ; num_AMs];
    time_arr = [time_arr; processing_time];
end
plot(num_AMs_arr, time_arr/60,'o-');

grid on
axis tight
hold off
xlabel("Number of Modules")
ylabel("time [min]")
title("Processing Time")
legend("Prev. Full", "New Full","New (9,5)","New (11,6)" ,'FontSize', 12)
