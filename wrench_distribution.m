num_down = 1;
L = [ones(1,5)*(lp+lg), ones(1,num_down)*(lp+lg+d)];
num_up = length(L) - num_down;
T_d = [];

N = length(x);
x = optimal(1:N*N_st);
u = optimal(N*N_st+1:end); 
x = reshape(x,N,N_st);%row vectors
u = reshape(u,N,N_u);%row vectors
t = linspace(0,1,N)'*t_f;

rows_per_matrix = 2; 
cols_per_matrix = 4*length(L); 
A_d = zeros(N * rows_per_matrix, N * cols_per_matrix);  % Preallocate full zero matrix


%% singularity
is_singularity = false(N,1);
for i = 1:N
    if(abs(x(i,2)-pi/2) <= 0.05 | abs(x(i,2)+pi/2) <= 0.05)
        is_singularity(i) = true;
        fprintf("\nat %f sec the robot is sigular \n",t(i))
    end
end
%%
for i= 1:N
    tau_ = get_tau_desired(x(i,:)',FDynamics(x(i,:)',u(i,:)'),L);
    if is_singularity(i)
         tau_=[0;0];%[0;tau_(2)];
    end    
    T_d = [T_d ; tau_];
    
    row_start = (i-1) * rows_per_matrix + 1;
    row_end = i * rows_per_matrix;
    col_start = (i-1) * cols_per_matrix + 1;
    col_end = i * cols_per_matrix;
 
    A_d(row_start:row_end, col_start:col_end) = get_At(x(i,:)',L);
end
Q = diag(ones(length(L)*4,1));
Q_d = kron(eye(N,N),Q);

%% QP
% Define the problem:
H = Q_d;  % Quadratic coefficient matrix (must be symmetric)
f = [];      % Linear coefficient vector
A = [];  % Inequality constraint matrix
b = [];              % Inequality constraint vector
Aeq = A_d;                    % equality constraints
beq = T_d;
%%%%%%% TODO %%%%%%
lb = [];%
ub = [];%[ ones(N*4*length(L) , 1) * 20];                  
thrust_bound = thrust_limit*Inf;
lb = repmat([ ones(4*num_up , 1) * 0 ; ones(4*num_down , 1) * -thrust_bound], N, 1) ;
ub = repmat([ ones(4*num_up , 1) * thrust_bound ; ones(4*num_down , 1) * 0], N, 1);                  
% Solve the QP problem using quadprog:
options = optimoptions('quadprog','Display','iter');  % Optional: to display iteration info
[U_d, fval, exitflag, output] = quadprog(H, f, A, b, Aeq, beq, lb, ub, []);

% Display results
%disp('Solution x:');
%disp(U_d);
disp('Objective function value at solution:');
disp(fval);

tau_error = Aeq*U_d - beq;
U_d_visual = reshape(U_d, [4*length(L),N])';
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