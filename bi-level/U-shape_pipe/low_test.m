clc
clear

r1 = [ 0  ,  0.8, 0];
r2 = [-0.8, -0.4, 0];
r3 = [ 0.8, -0.4, 0];

r = [r1, r2, r3];

mo = 2.0;
mi = 2.0;
g = 10.0;
f_max = 30;

Aeq = zeros(3,9);
beq = zeros(3,1);
Aeq(1:3,1:3) = eye(3);
Aeq(1:3,4:6) = eye(3);
Aeq(1:3,7:9) = eye(3);
beq(1:3) = (mo+3*mi)*g*[0;0;1];

x0 = ones(9, 1); % f1, f2, f3
y=[1,0,0];
%% gradient check
checkGradients(@(x)low_const(x), x0, Display="on", IsConstraint=true);
checkGradients(@(x)low_cost(x, y), x0, Display="on", IsConstraint=false);
%%
options = optimoptions('fmincon', 'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
  'display', 'iter-detailed', 'Algorithm','interior-point') ;  
[xsol, fval, exflag, output, multiplier, grad, hessian] = ...
  fmincon(@(x)low_cost(x, y), x0, [], [], Aeq, beq, [], [], @(x)low_const(x), options);

lambda = multiplier.ineqnonlin;
nu = multiplier.eqlin;

%%
skew(r1)*y'+2*lambda(1)*xsol(1:3)+nu
skew(r2)*y'+2*lambda(2)*xsol(4:6)+nu
skew(r3)*y'+2*lambda(3)*xsol(7:9)+nu
lambda(1) * (xsol(1:3)'*xsol(1:3) -900)
lambda(2) * (xsol(4:6)'*xsol(4:6) -900)
lambda(3) * (xsol(7:9)'*xsol(7:9) -900)
xsol(1:3)+xsol(4:6)+xsol(7:9)-[0;0;80];

function W = skew(w)
  W = zeros(3,3);
  W(1,2) = -w(3);
  W(1,3) = w(2);
  W(2,1) = w(3);
  W(2,3) = -w(1);
  W(3,1) = -w(2);
  W(3,2) = w(1);
end

function [obj, j] = low_cost(x, y_)  
r1 = [ 0  ,  0.8, 0];
r2 = [-0.8, -0.4, 0];
r3 = [ 0.8, -0.4, 0];
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