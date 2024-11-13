L = [lp+lg,lp+lg];
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

for i= 1:N
    T_d = [T_d ;  get_tau_desired(x(i,:)',FDynamics(x(i,:)',u(i,:)'),L)];
    
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
H = Q;  % Quadratic coefficient matrix (must be symmetric)
f = [];      % Linear coefficient vector
A = [];  % Inequality constraint matrix
b = [];              % Inequality constraint vector
lb = [];%
ub = [];%[ ones(N*4*length(L) , 1) * 20];                  
lb = [ ones(1*4*length(L) , 1) * 0];
ub = [ ones(1*4*length(L) , 1) * Inf];  
U_d = zeros(N, 4*length(L));
fval = zeros(N,1);
exitflag =zeros(N,1);

for i = 1:N
    Aeq = get_At(x(i,:)',L);                    % equality constraints
    beq = get_tau_desired(x(i,:)',FDynamics(x(i,:)',u(i,:)'),L);
    options = optimoptions('quadprog','Display','iter');  % Optional: to display iteration info
    [U_d(i,:), fval(i), exitflag(i)] = quadprog(H, f, A, b, Aeq, beq, lb, ub, []);
    
end

% Display results
disp('Objective function value at solution:');
disp(sum(fval));
%%
debug = A_d*reshape(U_d,[N*4*length(L),1]) - T_d;