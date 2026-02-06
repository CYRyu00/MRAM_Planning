clc
clear
close 

warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:illConditionedMatrix');

%% parameters
N = 50; % Number of points
y = fibonacci_sphere(N);

v1 = [ 0.02; 0.46; 0.0];
v2 = [ 0.01; 0.46; 0.01*sqrt(3)];
v3 = [-0.01; 0.46; 0.01*sqrt(3)];
v4 = [-0.02; 0.46; 0.0];
v5 = [-0.01; 0.46;-0.01*sqrt(3)];
v6 = [ 0.01; 0.46;-0.01*sqrt(3)];
v7 = [ 0.02;-0.46; 0.0];
v8 = [ 0.01;-0.46; 0.01*sqrt(3)];
v9 = [-0.01;-0.46; 0.01*sqrt(3)];
v10 = [-0.02;-0.46; 0.0];
v11 = [-0.01;-0.46;-0.01*sqrt(3)];
v12 = [ 0.01;-0.46;-0.01*sqrt(3)];

V = [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12];

sp1 = SupportFunction(V, 16);
R=[0,-1,0;1,0,0;0,0,1];
sp2 = SupportFunction(R*V, 16);
convex_shapes = {sp1, sp2, sp1};
c1 = [ 0.5;     1/6;0];
c2 = [ 0.0;-0.5+1/6;0];
c3 = [-0.5;     1/6;0];
centers = [c1,c2,c3];
object = ObjectGeometry(convex_shapes, centers);

% %% optimization variables
x1 = [1;0;0];
x2 = [0;-1;0];
x3 = [-1;0;0];
r1 = object.surface(1, x1);
r2 = object.surface(2, x2);
r3 = object.surface(3, x3);
z0 = [x1;x2;x3;r1;r2;r3];
% %%
% 
%% gradient check 
checkGradients(@(z)high_const(z, object), z0, Display="on", IsConstraint=true, Tolerance=1e-5);
checkGradients(@(z)high_cost(z, y, N), z0, Display="on", IsConstraint=false, Tolerance=1e-5);
%% 
options = optimoptions('fmincon', 'SpecifyConstraintGradient',true, 'SpecifyObjectiveGradient', true,...
  'display', 'iter-detailed', 'Algorithm','active-set');
tic
[zsol, fval, exflag, output, lambda, grad, hessian] = ...
  fmincon( @(z)high_cost(z, y, N), z0, [], [], [], [], [], [], @(z)high_const(z, object), options);
toc
%%
y = fibonacci_sphere(N);
[f,j] = high_cost(zsol,y,N);

function W = skew(w)
  W = zeros(3,3);
  W(1,2) = -w(3);
  W(1,3) = w(2);
  W(2,1) = w(3);
  W(2,3) = -w(1);
  W(3,1) = -w(2);
  W(3,2) = w(1);
end