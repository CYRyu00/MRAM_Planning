function plot_compare_results_ver2(x_opt1, u_opt1,dt,N,L,num_up1,x_opt2, u_opt2,num_up2)
%T = dt * N;
%% state
t = (0:1:N)*dt;
figure('Position',[50,250,600,500])
subplot(2,2,1)
plot(t, x_opt1(:,1),'.-')
hold on 
plot(t, x_opt2(:,1),'.-')
legend("Best","Worst")
ylabel('x [m]', 'Interpreter', 'latex')
title("States")
axis tight

subplot(2,2,2)
plot(t, x_opt1(:,2),'.-')
hold on 
plot(t, x_opt2(:,2),'.-')
legend("Best","Worst")
ylabel('$\theta$ [rad]', 'Interpreter', 'latex')

axis tight
subplot(2,2,3)
plot(t, x_opt1(:,3),'.-')
hold on 
plot(t, x_opt2(:,3),'.-')
legend("Best","Worst")
ylabel('$\dot{x}$ [m/s]', 'Interpreter', 'latex')
axis tight

subplot(2,2,4)
plot(t, x_opt1(:,4),'.-')
hold on 
plot(t, x_opt2(:,4),'.-')
legend("Best","Worst")
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex')


axis tight;

%% Best shape 
figure('Position',[50,50,800,400])
t = (1:1:N)*dt;
subplot(2,1,1)
hold on
legendEntries= [];
for i=1:length(num_up1)
    plot(t, u_opt1(:,2*i-1),'.-')
    legendEntries{i} = sprintf("AM %d", i); 
end
ylabel("f1 & f4 [N]")
legend(legendEntries);
title("Inputs : Best shape")
axis tight
hold off

subplot(2,1,2)
hold on
legendEntries= [];
for i=1:length(num_up1)
    plot(t, u_opt1(:,2*i),'.-')
    legendEntries{i} = sprintf("AM %d", i); 
end
ylabel("f2 & f3 [N]")
legend(legendEntries);
axis tight
hold off

%% 
taus1 = zeros(length(num_up1),N,2);

for i=1:length(num_up1)
    for j= 1:N
        tau_ = get_tau_distance(x_opt1(j,:)',u_opt1(j,2*i-1:2*i)',L(i));
        taus1(i,j,:)= tau_';
    end
end

t = (1:1:N)*dt;
sum_taus1 = sum(taus1, 1); 

%% Worst shape

figure('Position',[800,50,800,400])
t = (1:1:N)*dt;
subplot(2,1,1)
hold on
legendEntries= [];
for i=1:length(num_up2)
    plot(t, u_opt2(:,2*i-1),'.-')
    legendEntries{i} = sprintf("AM %d", i); 
end
ylabel("f1 & f4 [N]")
legend(legendEntries);
title("Inputs : Worst shape")
axis tight
hold off

subplot(2,1,2)
hold on
legendEntries= [];
for i=1:length(num_up2)
    plot(t, u_opt2(:,2*i),'.-')
    legendEntries{i} = sprintf("AM %d", i); 
end
ylabel("f2 & f3 [N]")
legend(legendEntries);
axis tight
hold off
%%
taus2 = zeros(length(num_up2),N,2);

for i=1:length(num_up2)
    for j= 1:N
        tau_ = get_tau_distance(x_opt1(j,:)',u_opt2(j,2*i-1:2*i)',L(i));
        taus2(i,j,:)= tau_';
    end
end

t = (1:1:N)*dt;
sum_taus2 = sum(taus2, 1); 

%% compare Best and Worst
figure('Position',[650,250,600,500])
subplot(2,1,1)
hold on
plot(t, sum_taus1(1,:,1),'.-');
plot(t, sum_taus2(1,:,1),'.-');
legend("Best","Worst");
title("Generalized Force: Best & Worst shape")
ylabel("tau1 [N]")
axis tight
hold off


subplot(2,1,2)
hold on
plot(t, sum_taus1(1,:,2),'.-');
plot(t, sum_taus2(1,:,2),'.-');
ylabel("tau2 [N·M]")
legend("Best","Worst");
axis tight
hold off
