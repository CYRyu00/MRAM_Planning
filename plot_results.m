N_st = 4;N_u=4;

%x_opt=[];
if ~isempty(x_opt)
    x=x_opt;
    u=u_opt;
    t = (0:1:N)*dt;

else
%shape 1
N = (length(optimal))/(N_st+N_u);
x = optimal(1:N*N_st);
u = optimal(N*N_st+1:end); 
x = reshape(x,N,N_st);%row vectors
u = reshape(u,N,N_u);%row vectors
t = (1:1:N)*dt;

end

%% state

figure(2)
subplot(2,2,1)
plot(t, x(:,1),'.-')
ylabel("x1")
axis tight

subplot(2,2,2)
plot(t, x(:,2),'.-')
ylabel("x2")

axis tight
subplot(2,2,3)
plot(t,x(:,3),'.-')
ylabel("x3")
axis tight

subplot(2,2,4)
plot(t, x(:,4),'.-')
ylabel("x4")
axis tight;

%% AM1
t = (0:1:N-1)*dt;
figure(3)
subplot(2,2,1)
plot(t, u(:,1),'.-')
ylabel("f1")
title("Inputs - AM1")
axis tight

subplot(2,2,2)
plot(t, u(:,2),'.-')
ylabel("f2")

axis tight
subplot(2,2,3)
plot(t,u(:,3),'.-')
ylabel("f3")

axis tight
subplot(2,2,4)
plot(t, u(:,4),'.-')
ylabel("f4")
axis tight;


%%
taus = zeros(N,2);
Fs = zeros(N,6);

for i= 1:N
    tau = get_tau(x(i,:)',u(i,:)');
    taus(i,:)= tau';
end
%%

figure(4);
subplot(2,1,1)
plot(t,taus(:,1),'.-');

title("tau1")
subplot(2,1,2)
plot(t,taus(:,2),'.-');

title("tau2")