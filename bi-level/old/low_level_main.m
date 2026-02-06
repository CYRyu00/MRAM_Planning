clc
clear

warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

mo = 2.0;
mi = 2.0;
g = 10.0;
f_max = 30;

r1 = [ 0  ;  0.8; 0];
r2 = [-0.8; -0.4; 0];
r3 = [ 0.8; -0.4; 0];
y0 = [1,0,0];
mu0 = [0,0,0];
z0 = [y0, mu0];

checkGradients(@(z)low_level_const(z), z0, Display="on", IsConstraint=true);
checkGradients(@(z)low_level_cost(z), z0, Display="on", IsConstraint=false);

%%

options = optimoptions('fmincon', 'SpecifyObjectiveGradient',true, 'SpecifyConstraintGradient',true, ...
  'display', 'iter-detailed', 'Algorithm','interior-point') ;

[zsol, fval, exflag, output, lambda, grad, hessian] = ...
  fmincon('low_level_cost', z0, [], [], [], [], [], [], 'low_level_const', options);

%%
zsol
fval