clc
clear
close 

warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:illConditionedMatrix');

%% parameters
N = 50; % Number of points
y = fibonacci_sphere(N);

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

sp = SupportFunction(V, 50);

%% optimization variables
x1 = [0;1;0];
x2 = [-1;0;0];
x3 = [1;0;0];
r1 = [ 0  ;  0.8; 0];
r2 = [-0.8; -0.4; 0];
r3 = [ 0.8; -0.4; 0];
z0 = [x1;x2;x3;r1;r2;r3];
%%

%% gradient check 
% checkGradients(@(z)high_const(z, sp), z0, Display="on", IsConstraint=true);
%clc
%[valid, err]=checkGradients(@(z)high_cost(z, y, N), z0, Display="on", IsConstraint=false);
%% 
options = optimoptions('fmincon', 'SpecifyConstraintGradient',true, 'SpecifyObjectiveGradient', true,...
  'display', 'iter-detailed', 'Algorithm','active-set');
tic
[zsol, fval, exflag, output, lambda, grad, hessian] = ...
  fmincon( @(z)high_cost(z, y, N), z0, [], [], [], [], [], [], @(z)high_const(z, sp), options);
toc

%%
% y = fibonacci_sphere(N);
% [f,j] = high_cost(zsol,y,N);

% function W = skew(w)
%   W = zeros(3,3);
%   W(1,2) = -w(3);
%   W(1,3) = w(2);
%   W(2,1) = w(3);
%   W(2,3) = -w(1);
%   W(3,1) = -w(2);
%   W(3,2) = w(1);
% end