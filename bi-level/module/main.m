clc
clear
close

warning('off', 'MATLAB:nearlySingularMatrix');

%% object
object = get_u_shape_pipe();
%% fibonacci sphere
N = 50; % Number of points
y = fibonacci_sphere(N);
%% init variables
n = 3;
x1 = [1;0;0];
x2 = [0;-1;0];
x3 = [-1;0;0];
r1 = object.surface(1, x1);
r2 = object.surface(2, x2);
r3 = object.surface(3, x3);
z0 = [x1;x2;x3;r1;r2;r3];
%% gradient check
checkGradients(@(z)high_const(z, object), z0, Display="on", IsConstraint=true, Tolerance=1e-5);
checkGradients(@(z)high_cost(z, y, N), z0, Display="on", IsConstraint=false, Tolerance=1e-5);
%% optimization
options = optimoptions('fmincon', 'SpecifyConstraintGradient',true, 'SpecifyObjectiveGradient', true,...
  'display', 'iter-detailed', 'Algorithm','active-set');
tic
[zsol, fval, exflag, output, lambda, grad, hessian] = ...
  fmincon( @(z)high_cost(z, y, N), z0, [], [], [], [], [], [], @(z)high_const(z, object), options);
toc
%% unpack
ox1 = zsol(1:3, 1); ox2 = zsol(4:6, 1); ox3 = zsol(7:9, 1);
or1 = zsol(10:12, 1); or2 = zsol(13:15, 1); or3 = zsol(16:18, 1);