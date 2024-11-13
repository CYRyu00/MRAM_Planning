function plot_multi_results(x_opt, u_opt,dt,N,L,num_up)
T = dt * N;
% state
t = (0:1:N)*dt;
figure
subplot(2,2,1)
plot(t, x_opt(:,1),'.-')
ylabel("x1")
title("State")
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
figure
t = (1:1:N)*dt;
subplot(2,1,1)
hold on
for i=1:length(num_up)
    plot(t, u_opt(:,2*i-1),'.-')
    legendEntries{i} = sprintf("AM %d", i); 
end
ylabel("f1 & f4")
legend(legendEntries);
title("Inputs")
axis tight
hold off

subplot(2,1,2)
hold on
for i=1:length(num_up)
    plot(t, u_opt(:,2*i),'.-')
    legendEntries{i} = sprintf("AM %d", i); 
end
ylabel("f2 & f3")
legend(legendEntries);
axis tight
hold off

%%
taus = zeros(length(num_up),N,2);
Fs = zeros(N,6);

for i=1:length(num_up)
    for j= 1:N
        tau_ = get_tau_distance(x_opt(j,:)',u_opt(j,2*i-1:2*i)',L(i));
        taus(i,j,:)= tau_';
    end
end

figure;
t = (1:1:N)*dt;
sum_taus = sum(taus, 1); 

subplot(2,1,1)
title("tau1")
hold on
for i=1:length(num_up)
    plot(t, taus(i,:,1),'.-');
    legendEntries{i} = sprintf("AM %d", i); 
end
legendEntries{i+1} = sprintf("Total");
plot(t, sum_taus(1,:,1),'.-');
legend(legendEntries);
axis tight
hold off


subplot(2,1,2)
title("tau2")
hold on
for i=1:length(num_up)
    plot(t, taus(i,:,2),'.-');
    legendEntries{i} = sprintf("AM %d", i); 
end
legendEntries{i+1} = sprintf("Total");
plot(t, sum_taus(1,:,2),'.-');
legend(legendEntries);
axis tight
hold off