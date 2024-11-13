global gridN
gridN = 100;
global t_f
t_f = 5;
global x_f
x_f = [3,0,0,0];

%Crazyflie System parameters
I = [16.571710 0.830806 0.718277;
     0.830806 16.655602 1.800197;
     0.718277 1.800197 29.261652]*10e-6;
thrust_limit = 0.15 ;
m0 = 0.028;
mu = 0.005964552;
r = 0.092*sqrt(2)/4;
d = 0.1;
lg = 0.05; % custom

m1 = 0.2; m2 = 0.1; lp = 0.2;
global params
params = [m1, m2, lp, lg, m0, I(2,2),mu,r];
%%

tic
% The initial parameter guess;  gridN states(4)  
% gridN inputs(4 thrust) 

X0 = [randn(gridN*8,1)];
% No linear inequality or equality constraints
A = [];
b = [];
Aeq = [];
Beq = [];
% lower bound and upper bound
% no constraints for states 
% thrust constraints unidirectional thrust
lb = [ ones(gridN*1 , 1) * -Inf ;  ones(gridN*4, 1) *  thrust_limit*0];
ub = [ ones(gridN*4, 1) *  Inf ;  ones(gridN*4, 1) *  Inf];

% Options for fmincon
options = optimoptions(@fmincon,'MaxFunctionEvaluations',1e6,...
                        'MaxIterations', 1000 ,...
                        'Display', 'iter', ...
                       'Algorithm', 'sqp');

% Solve for the best simulation time + control input
[optimal,fval,exitflag,output] = fmincon(@obj_func, X0, A, b, Aeq, Beq, lb, ub, ...
              @cons_func, options);
toc

fprintf("\n Time : %f", toc);
fprintf("\n Exitflag: %d, Cost: %4f\n",exitflag,fval)
