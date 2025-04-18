%close all

figure('Position',[300 500 500 400])
hold on
legend_entries = {};
max_iter_arr = [50,100,200,300,400,500];
num_AMs_arr = 5:1:12;
for MAX_ITER = max_iter_arr
    opt_val_arr=[];
    
    for num_AMs = num_AMs_arr
        filename = sprintf('result_9_5/max_iter_%d/%d_1_3.mat',MAX_ITER, num_AMs);
        load(filename);
        fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, processing_time)
        
        opt_val_arr = [opt_val_arr; optimal_value];
    end
    plot(num_AMs_arr, opt_val_arr,'o-');
    legend_entries{end+1} = sprintf('%d', MAX_ITER);
end
legend(legend_entries)
xlabel("Number of Modules")
ylabel("J*")
title("object function value")
grid on
%
figure('Position',[800 500 500 400])
hold on
legend_entries = {};
for MAX_ITER = max_iter_arr
    
    time_arr=[];
    
    for num_AMs = num_AMs_arr
        filename = sprintf('result_9_5/max_iter_%d/%d_1_3.mat',MAX_ITER, num_AMs);
        load(filename);
        fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, processing_time)
        
        time_arr = [time_arr; processing_time];
    end
    plot(num_AMs_arr, time_arr,'o-');
    legend_entries{end+1} = sprintf('%d', MAX_ITER);
end
legend(legend_entries)
xlabel("Number of Modules")
ylabel("[sec]")
title("Processing time")
grid on
%%
%close all

figure('Position',[300 500 500 400])
hold on
legend_entries = {};
gamma_arr = 0.2:0.1:0.5;
num_AMs_arr = 5:1:12;
for GAMMA = gamma_arr
    opt_val_arr=[];
    
    for num_AMs = num_AMs_arr
        filename = sprintf('result_9_5/max_iter_400/gamma_%.2f/%d_1_3.mat', GAMMA, num_AMs);
        load(filename);
        fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, processing_time)
        
        opt_val_arr = [opt_val_arr; optimal_value];
    end
    plot(num_AMs_arr, opt_val_arr,'o-');
    legend_entries{end+1} = sprintf('%.2f',GAMMA);
end
legend(legend_entries)
xlabel("Number of Modules")
ylabel("J*")
title("object function value")
grid on
%
figure('Position',[800 500 500 400])
hold on
legend_entries = {};
for GAMMA = gamma_arr
    
    time_arr=[];
    
    for num_AMs = num_AMs_arr
        filename = sprintf('result_9_5/max_iter_400/gamma_%.2f/%d_1_3.mat', GAMMA, num_AMs);
        load(filename);
        fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, processing_time)
        
        time_arr = [time_arr; processing_time];
    end
    plot(num_AMs_arr, time_arr,'o-');
    legend_entries{end+1} = sprintf('%.2f', GAMMA);
end
legend(legend_entries)
xlabel("Number of Modules")
ylabel("[sec]")
title("Processing time")
grid on
%%
%close all

figure('Position',[300 500 500 400])
hold on
legend_entries = {};
eps_arr = 0.05:0.05:0.25;
num_AMs_arr = 5:1:12;
for EPS = eps_arr
    opt_val_arr=[];
    
    for num_AMs = num_AMs_arr
        filename = sprintf('result_9_5/max_iter_400/gamma_0.20/eps_%.2f/%d_1_3.mat', EPS, num_AMs);
        load(filename);
        fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, processing_time)
        
        opt_val_arr = [opt_val_arr; optimal_value];
    end
    plot(num_AMs_arr, opt_val_arr,'o-');
    legend_entries{end+1} = sprintf('%.2f',EPS);
end
legend(legend_entries)
xlabel("Number of Modules")
ylabel("J*")
title("object function value")
grid on
%
figure('Position',[800 500 500 400])
hold on
legend_entries = {};
for EPS = eps_arr
    
    time_arr=[];
    
    for num_AMs = num_AMs_arr
        filename = sprintf('result_9_5/max_iter_400/gamma_0.20/eps_%.2f/%d_1_3.mat', EPS, num_AMs);
        load(filename);
        fprintf("num AMs: %d, J*: %f, processing time: %f\n", num_AMs, optimal_value, processing_time)
        
        time_arr = [time_arr; processing_time];
    end
    plot(num_AMs_arr, time_arr,'o-');
    legend_entries{end+1} = sprintf('%.2f', EPS);
end
legend(legend_entries)
xlabel("Number of Modules")
ylabel("[sec]")
title("Processing time")
grid on