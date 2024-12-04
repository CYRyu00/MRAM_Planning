function plot_compare_results(x_opt1, u_opt1,dt,N,L,num_up1,x_opt2, u_opt2,num_up2)
%T = dt * N;
%% state
t = (0:1:N)*dt;
figure('Position',[50,250,600,500])
subplot(2,2,1)
plot(t, x_opt1(:,1),'.-')
hold on 
plot(t, x_opt2(:,1),'.-')
legend("1st","2nd")
ylabel("x [m/s]")
title("States")
axis tight

subplot(2,2,2)
plot(t, x_opt1(:,2),'.-')
hold on 
plot(t, x_opt2(:,2),'.-')
legend("1st","2nd")
ylabel("theta [rad]")

axis tight
subplot(2,2,3)
plot(t, x_opt1(:,3),'.-')
hold on 
plot(t, x_opt2(:,3),'.-')
legend("1st","2nd")
ylabel("x\_dot [m/s^2]")
axis tight

subplot(2,2,4)
plot(t, x_opt1(:,4),'.-')
hold on 
plot(t, x_opt2(:,4),'.-')
legend("1st","2nd")
ylabel("theta\_dot [rad/s]")
axis tight;

%% 1st shape 
figure('Position',[50,50,800,400])
t = (1:1:N)*dt;
subplot(2,2,1)
hold on
legendEntries= [];
for i=1:length(num_up1)
    plot(t, u_opt1(:,2*i-1),'.-')
    legendEntries{i} = sprintf("AM %d", i); 
end
ylabel("f1 & f4 [N]")
legend(legendEntries);
title("Inputs : 1st shape")
axis tight
hold off

subplot(2,2,2)
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

subplot(2,2,3)
hold on
legendEntries= [];
for i=1:length(num_up1)
    plot(t, taus1(i,:,1),'.-');
    legendEntries{i} = sprintf("AM %d", i); 
end
legendEntries{i+1} = sprintf("Total");
plot(t, sum_taus1(1,:,1),'.-');
legend(legendEntries);
title("Generalized Force : 1st shape")
ylabel("tau1 [N]")
axis tight
hold off


subplot(2,2,4)
hold on
legendEntries= [];
for i=1:length(num_up1)
    plot(t, taus1(i,:,2),'.-');
    legendEntries{i} = sprintf("AM %d", i); 
end
legendEntries{i+1} = sprintf("Total");
plot(t, sum_taus1(1,:,2),'.-');
legend(legendEntries);
ylabel("tau2 [N·M]")
axis tight
hold off

%% 2nd shape

figure('Position',[800,50,800,400])
t = (1:1:N)*dt;
subplot(2,2,1)
hold on
legendEntries= [];
for i=1:length(num_up2)
    plot(t, u_opt2(:,2*i-1),'.-')
    legendEntries{i} = sprintf("AM %d", i); 
end
ylabel("f1 & f4 [N]")
legend(legendEntries);
title("Inputs : 2nd shape")
axis tight
hold off

subplot(2,2,2)
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

subplot(2,2,3)
hold on
legendEntries= [];
for i=1:length(num_up2)
    plot(t, taus2(i,:,1),'.-');
    legendEntries{i} = sprintf("AM %d", i); 
end
legendEntries{i+1} = sprintf("Total");
plot(t, sum_taus2(1,:,1),'.-');
legend(legendEntries);
title("Generalized Force : 2nd shape")
ylabel("tau1 [N]")
axis tight
hold off

subplot(2,2,4)
hold on
legendEntries= [];
for i=1:length(num_up2)
    plot(t, taus2(i,:,2),'.-');
    legendEntries{i} = sprintf("AM %d", i); 
end
legendEntries{i+1} = sprintf("Total");
plot(t, sum_taus2(1,:,2),'.-');
legend(legendEntries);
ylabel("tau2 [N·M]")
axis tight
hold off

%% compare 1st and 2nd
figure('Position',[650,250,600,500])
subplot(2,1,1)
hold on
plot(t, sum_taus1(1,:,1),'.-');
plot(t, sum_taus2(1,:,1),'.-');
legend("1st total","2nd total");
title("Generalized Force: 1st & 2nd shape")
ylabel("tau1 [N]")
axis tight
hold off


subplot(2,1,2)
hold on
plot(t, sum_taus1(1,:,2),'.-');
plot(t, sum_taus2(1,:,2),'.-');
ylabel("tau2 [N·M]")
legend("1st total","2nd total");
axis tight
hold off
