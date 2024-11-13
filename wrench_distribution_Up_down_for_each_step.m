num_down = 1;
L = [ones(1,5)*(lp+lg), ones(1,num_down)*(lp+lg)];
num_up = length(L) - num_down;
T_d = [];

N = length(x);
x = optimal(1:N*N_st);
u = optimal(N*N_st+1:end); 
x = reshape(x,N,N_st);%row vectors
u = reshape(u,N,N_u);%row vectors
t = linspace(0,1,N)'*t_f;

for i= 1:N
    T_d = [T_d ;  get_tau_desired(x(i,:)',FDynamics(x(i,:)',u(i,:)'),L)];
end
Q = diag(ones(length(L)*4,1));

%% singularity
is_singularity = false(N,1);
for i = 1:N
    if(abs(x(i,2)-pi/2) <= 0.02)
        is_singularity(i) = true;
        fprintf("\nat %f sec the robot is sigular \n",t(i))
    end
end
%% QP
% Define the problem:
H = Q;       % Quadratic coefficient matrix (must be symmetric)
f = [];      % Linear coefficient vector
A = [];      % Inequality constraint matrix
b = [];              % Inequality constraint vector
lb = [];
ub = [];
thrust_bound = thrust_limit*Inf;
lb = [ ones(4*num_up , 1) * 0 ; ones(4*num_down , 1) * - thrust_bound];
ub = [ ones(4*num_up , 1) * thrust_bound ; ones(4*num_down , 1) * 0];  

U_d = zeros(N, 4*length(L));
fval = zeros(N,1);
exitflag =zeros(N,1);

for i = 1:N
    Aeq = get_At(x(i,:)',L);                  
    beq = get_tau_desired(x(i,:)',FDynamics(x(i,:)',u(i,:)'),L);
    if is_singularity(i)
        beq = beq(2);
        Aeq = Aeq(2,:);
    end
    options = optimoptions('quadprog','Display','iter');  % Optional: to display iteration info
    [U_d(i,:), fval(i), exitflag(i)] = quadprog(H, f, A, b, Aeq, beq, lb, ub, []);
    
end
U_d_visual = U_d;
T_d_visual = reshape(T_d ,[2,N])';
%% INputs
figure(3)
subplot(2,2,1)
plot(t, u(:,1),'.-')
hold on
for i=1:length(L)
    plot(t, U_d_visual(:,1 + (i-1)*4),'.-')
end 
hold off
legend("1AM\_wo\_cons","AM1","AM2","AM3","AM4");
ylabel("f1")
title("Inputs - AM1")
axis tight

subplot(2,2,2)
plot(t, u(:,2),'.-')
hold on
for i=1:length(L)
    plot(t, U_d_visual(:,2 + (i-1)*4),'.-')
end 
hold off
legend("1AM\_wo\_cons","AM1","AM2","AM3","AM4");
ylabel("f2")

axis tight
subplot(2,2,3)
plot(t,u(:,3),'.-')
hold on
for i=1:length(L)
    plot(t, U_d_visual(:,3 + (i-1)*4),'.-')
end 
hold off
legend("1AM\_wo\_cons","AM1","AM2","AM3","AM4");ylabel("f3")

axis tight
subplot(2,2,4)
plot(t, u(:,4),'.-')
ylabel("f4")
hold on
for i=1:length(L)
    plot(t, U_d_visual(:,4 + (i-1)*4),'.-')
end 
hold off
legend("1AM\_wo\_cons","AM1","AM2","AM3","AM4");
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
plot(t,taus(:,1),'.-',t,T_d_visual(:,1),'.-');
legend("1AM\_wo\_cons", "Multi\_AM")

title("tau1")
subplot(2,1,2)
plot(t,taus(:,2),t,T_d_visual(:,2),'.-');
legend("1AM\_wo\_cons", "Multi\_AM")
title("tau2")