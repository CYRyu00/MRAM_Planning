clc
clear
close 

%% parameters
N = 20; % Number of points
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

sp = SupportFunction(V, 20);

%% optimization variables
x1 = [0;1;0];
x2 = [-1;0;0];
x3 = [1;0;0];
r1 = [ 0  ;  0.8; 0];
r2 = [-0.8; -0.5; 0];
r3 = [ 0.8; -0.4; 0];
z0 = [x1;x2;x3;r1;r2;r3];
%%

%% gradient check 
% checkGradients(@(z)high_const(z, sp), z0, Display="on", IsConstraint=true);
clc
[f,j] = tmp_cost(z0);
%%
clc
[valid, err]=checkGradients(@(z)tmp_cost(z), z0, Display="on", IsConstraint=false);

%% 
% options = optimoptions('fmincon', 'SpecifyConstraintGradient',true, 'SpecifyObjectiveGradient',true, ...
%   'display', 'iter-detailed', 'Algorithm','active-set') ;
% 
% [zsol, fval, exflag, output, lambda, grad, hessian] = ...
%   fmincon( @(z)tmp_cost(z), z0, [], [], [], [], [], [], @(z)high_const(z, sp), options);

%%
% r1 = zsol(10:12, 1)
% r2 = zsol(13:15, 1)
% r3 = zsol(16:18, 1)


function [obj, j] = tmp_cost(z)
  y = [sqrt(0.5), sqrt(0.5),0];
  r1 = z(10:12, 1);
  r2 = z(13:15, 1);
  r3 = z(16:18, 1);
  yk = y;

  Aeq = zeros(3,9);
  beq = zeros(3,1);
  Aeq(1:3,1:3) = eye(3);
  Aeq(1:3,4:6) = eye(3);
  Aeq(1:3,7:9) = eye(3);
  beq(1:3) = [0;0;80];
  
  x0 = ones(9, 1); % f1, f2, f3

  options = optimoptions('fmincon', 'SpecifyObjectiveGradient',true,...
    'SpecifyConstraintGradient',true,...
    'Display', 'none', 'Algorithm','interior-point') ;
  [xsol, fval, ~, ~, multiplier] = ...
    fmincon(@(x)low_cost(x, yk), x0, [], [], Aeq, beq, [], [],...
    @(x)low_const(x), options);
  obj = fval;
  lambda = multiplier.ineqnonlin;


  if nargout >1
    j = zeros(18, 1);
    Q = zeros(15, 15);
    b = zeros(15, 9);

    lambdak = lambda;
    xk = xsol;

    Q(1:3,1:3) = 2 * lambdak(1)*eye(3);
    Q(4:6,4:6) = 2 * lambdak(2)*eye(3);
    Q(7:9,7:9) = 2 * lambdak(3)*eye(3);

    Q(1:3,10) = 2 * xk(1:3);
    Q(4:6,11) = 2 * xk(4:6);
    Q(7:9,12) = 2 * xk(7:9);
    Q(10,1:3) = 2 * lambdak(1)*xk(1:3);
    Q(11,4:6) = 2 * lambdak(2)*xk(4:6);
    Q(12,7:9) = 2 * lambdak(3)*xk(7:9);
    Q(10, 10) = xk(1:3)'*xk(1:3)-900;
    Q(11, 11) = xk(4:6)'*xk(4:6)-900;
    Q(12, 12) = xk(7:9)'*xk(7:9)-900;

    Q(1:3,13:15) = eye(3);
    Q(4:6,13:15) = eye(3);
    Q(7:9,13:15) = eye(3);
    Q(13:15,1:3) = eye(3);
    Q(13:15,4:6) = eye(3);
    Q(13:15,7:9) = eye(3);

    b(1:3,1:3) = skew(yk);
    b(4:6,4:6) = skew(yk);
    b(7:9,7:9) = skew(yk);

    dxdr = Q\b;
    dfdr = dxdr(1:9,1:9);
    depsilondr = zeros(1,9);
    depsilondr(1:3) = -yk*(skew(xk(1:3)-[0;0;20])) ;
    depsilondr(4:6) = -yk*(skew(xk(4:6)-[0;0;20])) ;
    depsilondr(7:9) = -yk*(skew(xk(7:9)-[0;0;20])) ;
    depsilondr = depsilondr + yk*skew(r1) * dfdr(1:3,1:9);
    depsilondr = depsilondr + yk*skew(r2) * dfdr(4:6,1:9);
    depsilondr = depsilondr + yk*skew(r3) * dfdr(7:9,1:9);
    j(10:18) = depsilondr';
    j = -j;
  end

  function [obj, j] = low_cost(x, y_)  
    tau1 = skew(r1) * (x(1:3) - [0;0;20]);
    tau2 = skew(r2) * (x(4:6) - [0;0;20]);
    tau3 = skew(r3) * (x(7:9) - [0;0;20]);
  
    tau = tau1+tau2+tau3;
    obj = -y_*tau;
    if nargout >1
      j = zeros(9, 1);
      j(1:3) = -y_*skew(r1);
      j(4:6) = -y_*skew(r2);
      j(7:9) = -y_*skew(r3);
    end
  end
  
  function [g,h,gg,gh] = low_const(x)
    g(1) = x(1:3)'*x(1:3)-900;
    g(2) = x(4:6)'*x(4:6)-900;
    g(3) = x(7:9)'*x(7:9)-900;
    h = [];
  
    if nargout > 2
      gg = zeros(9,3);
      gg(1:3, 1) = 2*x(1:3);
      gg(4:6, 2) = 2*x(4:6);
      gg(7:9, 3) = 2*x(7:9);
  
      gh = [];
    end
  end
end

function W = skew(w)
  W = zeros(3,3);
  W(1,2) = -w(3);
  W(1,3) = w(2);
  W(2,1) = w(3);
  W(2,3) = -w(1);
  W(3,1) = -w(2);
  W(3,2) = w(1);
end
