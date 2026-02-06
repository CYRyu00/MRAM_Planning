% clc
% clear

% r1 = [ 0  ,  0.8, 0];
% r2 = [-0.8, -0.4, 0];
% r3 = [ 0.8, -0.4, 0];
% zsol = [0,0,0,0,0,0,0,0,0,r1,r2,r3]
u0 = [sqrt(0.5),sqrt(0.5),0];
options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true, 'display', 'iter-detailed','Algorithm','interior-point');
% options = optimoptions('fmincon', 'display', 'none', 'Algorithm','interior-point') ;  
[dir_sol, radius] = fmincon(@(u)high_level_cost(u, zsol), u0, [], [], [], [], [], [], @(u)high_level_constraint(u), options);

dir_sol
radius

function [obj, j] = high_level_cost(u, zsol)
  r1 = zsol(10:12);
  r2 = zsol(13:15);
  r3 = zsol(16:18);
  mo = 2.0;
  mi = 2.0;
  g = 10.0;

  Aeq = zeros(6,10);
  beq = zeros(6,1);
  Aeq(1:3,1) = u';
  Aeq(1:3,2:4) = -skew(r1);
  Aeq(1:3,5:7) = -skew(r2);
  Aeq(1:3,8:10) = -skew(r3);
  Aeq(4:6,2:4) = eye(3);
  Aeq(4:6,5:7) = eye(3);
  Aeq(4:6,8:10) = eye(3);
  beq(1:3) = -mi * g * (skew(r1+r2+r3)*[0;0;1]);
  beq(4:6) = (mo+3*mi)*g*[0;0;1];
  x0 = zeros(10,1);
  options = optimoptions('fmincon','display','none', 'Algorithm','interior-point') ;  
  [sol, ~, ~, ~, labmda] = fmincon(@(x)low_level_cost(x), x0, [], [], Aeq, beq, [], [], @(x)low_level_constraint(x), options);
  
  lambda_g = labmda.ineqnonlin;

  nu_1 = labmda.eqlin(1:3);
  nu_2 = labmda.eqlin(4:6);
  obj = sol(1);
  if nargout >1
    Q = zeros(19, 19);
    b = zeros(19, 3);

    Q(1, 14:16) = u;
    Q(14:16, 1) = u';

    Q(2:4,2:4) = lambda_g(1) * eye(3);
    Q(5:7,5:7) = lambda_g(2) * eye(3);
    Q(8:10,8:10) = lambda_g(3) * eye(3);

    Q(2:4,11) = sol(2:4)';
    Q(5:7,12) = sol(5:7)';
    Q(8:10,13) = sol(8:10)';
    Q(11,2:4) = lambda_g(1)*sol(2:4);
    Q(12,5:7) = lambda_g(2)*sol(5:7);
    Q(13,8:10) = lambda_g(3)*sol(8:10);

    Q(11,11) = sol(2:4)'*sol(2:4)-900;
    Q(12,12) = sol(5:7)'*sol(5:7)-900;
    Q(13,13) = sol(8:10)'*sol(8:10)-900;

    Q(2:4,14:16) = skew(r1);
    Q(5:7,14:16) = skew(r2);
    Q(8:10,14:16) = skew(r3);
    Q(14:16,2:4) = -skew(r1);
    Q(14:16,5:7) = -skew(r2);
    Q(14:16,8:10) = -skew(r3);

    Q(2:4,17:19) = eye(3);
    Q(5:7,17:19) = eye(3);
    Q(8:10,17:19) = eye(3);
    Q(17:19,2:4) = eye(3);
    Q(17:19,5:7) = eye(3);
    Q(17:19,8:10) = eye(3);
    
    b(1,1:3) = -nu_1;
    b(14:16,1:3) = -sol(1)*eye(3);

    dzdu = Q\b;
    j = dzdu(1,1:3);
  end
end
function [g,h,gg,gh] = high_level_constraint(u)
  g =[];
  h = u*u'-1;
  if nargout >2
    gg = [];
    gh = 2*u';
  end
end

function f = low_level_cost(x)
 f = -x(1);
end
function [g,h,gg,gh] = low_level_constraint(x)
  g(1) = x(2:4)'*x(2:4)-900;
  g(2) = x(5:7)'*x(5:7)-900;
  g(3) = x(8:10)'*x(8:10)-900;
  h = [];

  if nargout >2
    gg = zeros(3,10);
    gg(1, 2:4) = x(2:4);
    gg(2, 5:7) = x(5:7);
    gg(3, 8:10) = x(8:10);
    gh = [];
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