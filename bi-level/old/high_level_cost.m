function [obj, j] = high_level_cost(z)
  u0 = [1,0,0];
  middle_options = optimoptions('fmincon', 'SpecifyObjectiveGradient',true, 'SpecifyConstraintGradient',true,...
    'display', 'none', 'Algorithm','interior-point') ;
  
  low_sol = []; % alpha, f1, f2, f3
  low_lambda = []; %

  [middle_sol, ~, ~, ~, middle_lambda, middle_grad, middle_hessian] = ...
    fmincon(@(u)middle_level_cost(u, z), u0, [], [], [], [], [], [], 'middle_level_const', middle_options);
  u0 = middle_sol;
  obj = -low_sol(1);

  if nargout >1    
    nu_2 = low_lambda.eqlin(4:6);
    lambda_g = low_lambda.ineqnonlin;
    u = middle_sol;
    lambda_m = middle_lambda.ineqnonlin;
    alpha = low_sol(1);
    r1 = z(10:12, 1);
    r2 = z(13:15, 1);
    r3 = z(16:18, 1);

    Q2 = zeros(23, 23);
    b2 = zeros(23, 9);

    Q2(1, 17:19) = u;
    Q2(17:19, 1) = u';
    Q2(1, 20:22) = nu_2;

    Q2(2:4,2:4) = 2*lambda_g(1) * eye(3);
    Q2(5:7,5:7) = 2*lambda_g(2) * eye(3);
    Q2(8:10,8:10) = 2*lambda_g(3) * eye(3);

    Q2(2:4,11) = 2*low_sol(2:4)';
    Q2(5:7,12) = 2*low_sol(5:7)';
    Q2(8:10,13) = 2*low_sol(8:10)';
    Q2(11,2:4) = 2*lambda_g(1)*low_sol(2:4);
    Q2(12,5:7) = 2*lambda_g(2)*low_sol(5:7);
    Q2(13,8:10) = 2*lambda_g(3)*low_sol(8:10);

    Q2(11,11) = low_sol(2:4)'*low_sol(2:4)-900;
    Q2(12,12) = low_sol(5:7)'*low_sol(5:7)-900;
    Q2(13,13) = low_sol(8:10)'*low_sol(8:10)-900;

    Q2(2:4,14:16) = eye(3);
    Q2(5:7,14:16) = eye(3);
    Q2(8:10,14:16) = eye(3);
    Q2(14:16,2:4) = eye(3);
    Q2(14:16,5:7) = eye(3);
    Q2(14:16,8:10) = eye(3);

    Q2(2:4,17:19) = skew(r1);
    Q2(5:7,17:19) = skew(r2);
    Q2(8:10,17:19) = skew(r3);
    Q2(17:19,2:4) = -skew(r1);
    Q2(17:19,5:7) = -skew(r2);
    Q2(17:19,8:10) = -skew(r3);
    Q2(17:19, 20:22) = alpha * eye(3);
    
    Q2(20:22, 1) = middle_grad;
    Q2(20:22, 20:22) = 2*lambda_m * eye(3)+middle_hessian;
    Q2(20:22, 23) = 2*u';
    Q2(23, 20:22) = 2*lambda_m*u;

    b2(2:4, 1:3) = -skew(nu_2);
    b2(5:7, 4:6) = -skew(nu_2);
    b2(8:10, 7:9) = -skew(nu_2);

    b2(17:19, 1:3) = skew(20*[0;0;1]-low_sol(2:4));
    b2(17:19, 4:6) = skew(20*[0;0;1]-low_sol(5:7));
    b2(17:19, 7:9) = skew(20*[0;0;1]-low_sol(8:10));

  
    dzdr = Q2\b2;

    dalphadr = dzdr(1,1:9);
    j = zeros(18, 1);
    j(10:18) = -dalphadr;
    j = j';
  end
  
  function [obj, j] = middle_level_cost(u, z)
    r1 = z(10:12, 1);
    r2 = z(13:15, 1);
    r3 = z(16:18, 1);
    
    mo = 2.0;
    mi = 2.0;
    g = 10.0;
    f_max = 30;
  
    Aeq = zeros(6,10);
    beq = zeros(6,1);
    Aeq(1:3,2:4) = eye(3);
    Aeq(1:3,5:7) = eye(3);
    Aeq(1:3,8:10) = eye(3);
    Aeq(4:6,1) = u';
    Aeq(4:6,2:4) = -skew(r1);
    Aeq(4:6,5:7) = -skew(r2);
    Aeq(4:6,8:10) = -skew(r3);
  
    beq(1:3) = (mo+3*mi)*g*[0;0;1];
    beq(4:6) = -mi * g * (skew(r1+r2+r3)*[0;0;1]);
  
    x0 = zeros(10,1);
    low_options = optimoptions('fmincon','display','none', 'Algorithm','interior-point', 'SpecifyObjectiveGradient',true, 'SpecifyConstraintGradient',true) ;  
    [low_sol, ~, ~, ~, low_lambda] = fmincon('low_level_cost', x0, [], [], Aeq, beq, [], [], 'low_level_const', low_options);
    obj = low_sol(1);

    nu_2 = low_lambda.eqlin(4:6);
    lambda_g = low_lambda.ineqnonlin;
    
    if nargout >1
      Q = zeros(19, 19);
      b = zeros(19, 3);

      Q(1, 17:19) = u';
      Q(17:19, 1) = u;
  
      Q(2:4,2:4) = 2*lambda_g(1) * eye(3);
      Q(5:7,5:7) = 2*lambda_g(2) * eye(3);
      Q(8:10,8:10) = 2*lambda_g(3) * eye(3);
  
      Q(2:4,11) = 2*low_sol(2:4)';
      Q(5:7,12) = 2*low_sol(5:7)';
      Q(8:10,13) = 2*low_sol(8:10)';
      Q(11,2:4) = 2*lambda_g(1)*low_sol(2:4);
      Q(12,5:7) = 2*lambda_g(2)*low_sol(5:7);
      Q(13,8:10) = 2*lambda_g(3)*low_sol(8:10);
  
      Q(11,11) = low_sol(2:4)'*low_sol(2:4)-900;
      Q(12,12) = low_sol(5:7)'*low_sol(5:7)-900;
      Q(13,13) = low_sol(8:10)'*low_sol(8:10)-900;
  
      Q(2:4,14:16) = eye(3);
      Q(5:7,14:16) = eye(3);
      Q(8:10,14:16) = eye(3);
      Q(14:16,2:4) = eye(3);
      Q(14:16,5:7) = eye(3);
      Q(14:16,8:10) = eye(3);
  
      Q(2:4,17:19) = skew(r1);
      Q(5:7,17:19) = skew(r2);
      Q(8:10,17:19) = skew(r3);
      Q(17:19,2:4) = -skew(r1);
      Q(17:19,5:7) = -skew(r2);
      Q(17:19,8:10) = -skew(r3);
  
      b(1,1:3) = -nu_2;
      b(17:19,1:3) = -low_sol(1)*eye(3);
  
      dzdu = Q\b;
      j = dzdu(1,1:3);
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