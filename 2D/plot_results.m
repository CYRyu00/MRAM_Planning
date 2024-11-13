function plot_results(x_opt, u_opt,dt,N)
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
plot(t, u_opt(:,1),'.-')
ylabel("f1 & f4")
title("Inputs")
axis tight

subplot(2,1,2)
plot(t, u_opt(:,2),'.-')
ylabel("f2 & f3")
%%
taus = zeros(N,2);
Fs = zeros(N,6);

for i= 1:N
    tau_ = get_tau(x_opt(i,:)',u_opt(i,:)');
    taus(i,:)= tau_';
end

figure;
subplot(2,1,1)
plot(t,taus(:,1),'.-');

title("tau1")
subplot(2,1,2)
plot(t,taus(:,2),'.-');

title("tau2")
