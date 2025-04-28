global num_down
num_down = 0;
global L
L = [ones(1,2)*(lp+lg), ones(1,num_down)*(lp+lg)];
num_up = length(L) - num_down;
global gridN

%% NLP
tic
% The initial parameter guess;  gridN states(4)  
% gridN inputs(4 thrust) 
N_st = 4;N_u=4;
N_ = (length(optimal))/(N_st+N_u);
x = optimal(1:N_*N_st);
u = optimal(N_*N_st+1:end); 
u = [ repmat(u/num_up,num_up,1);zeros(length(u)*num_down,1)];
X0 = [x;u];

% No linear inequality or equality constraints
A = [];
b = [];
Aeq = [];
Beq = [];
% lower bound and upper bound
% no constraints for states 
% thrust constraints unidirectional thrust
lb = [ ones(gridN*4 , 1) * -Inf ;  ... %states
       ones(gridN*4* num_up,1) *  thrust_limit*0; ... % upside AMs
       ones(gridN*4* num_down,1) *  thrust_limit*-Inf]; % downside AMs
ub = [ ones(gridN*4, 1) *  Inf ;  ... 
       ones(gridN*4* num_up,1) *  thrust_limit*Inf; ...
       ones(gridN*4* num_down,1) *  thrust_limit*0];

% Options for fmincon
options = optimoptions(@fmincon,'MaxFunctionEvaluations',1e6,...
                        'MaxIterations', 1000 ,...
                        'Display', 'iter', ...
                       'Algorithm', 'sqp');

% Solve for the best simulation time + control input
[new_optimal,fval,exitflag,output] = fmincon(@obj_func, X0, A, b, Aeq, Beq, lb, ub, ...
              @cons_func, options);
toc

fprintf("\n Time : %f", toc);
fprintf("\n Exitflag: %d, Cost: %4f\n",exitflag,fval)
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