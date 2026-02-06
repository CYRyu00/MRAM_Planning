clc
clear

warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

mo = 2.0;
mi = 2.0;
g = 10.0;
f_max = 30;

v1 = [ 0.5; 0.5; 0.02];
v2 = [ 0.5;-0.5; 0.02];
v3 = [ 0.0;-0.5; 0.02];
v4 = [-0.5;-0.5; 0.02];
v5 = [-0.5; 0.5; 0.02];
v6 = [ 0.0; 0.5; 0.02];
v7 = [ 0.5; 0.5;-0.02];
v8 = [ 0.5;-0.5;-0.02];
v9 = [ 0.0;-0.5;-0.02];
v10 = [-0.5;-0.5;-0.02];
v11 = [-0.5; 0.5;-0.02];
v12 = [ 0.0; 0.5;-0.02];

V = [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12];
global sp;
sp = SupportFunction(V, 10);

x1 = [0;1;0];
x2 = [-1;0;0];
x3 = [1;0;0];
r1 = [ 0  ;  0.8; 0];
r2 = [-0.8; -0.4; 0];
r3 = [ 0.8; -0.4; 0];
z0 = [x1;x2;x3;r1;r2;r3];

checkGradients(@(z)high_level_const(z), z0, Display="on", IsConstraint=true);
checkGradients(@(z)high_level_cost(z), z0, Display="on", IsConstraint=false);

%%

options = optimoptions('fmincon', 'SpecifyObjectiveGradient',true, 'SpecifyConstraintGradient',true, ...
  'display', 'iter-detailed', 'Algorithm','active-set') ;

[zsol, fval, exflag, output, lambda, grad, hessian] = ...
  fmincon('high_level_cost', z0, [], [], [], [], [], [], 'high_level_const', options);

%%