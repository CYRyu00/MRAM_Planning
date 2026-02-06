function [obj, j] = high_cost(z, y, N)
  r1 = z(10:12, 1);
  r2 = z(13:15, 1);
  r3 = z(16:18, 1);
  expsum = 0;
  beta = 5.0;
  w = zeros(N, 1);
  lambda = zeros(N, 3);
  nu = zeros(N, 6);
  x = zeros(N, 10);
  parfor k = 1:N    
    yk = y(k, :);

    Aeq = zeros(6,9);
    beq = zeros(6,1);
    Aeq(1:3,2:4)  = eye(3);
    Aeq(1:3,5:7)  = eye(3);
    Aeq(1:3,8:10) = eye(3);
    Aeq(4:6,1)  = yk';
    Aeq(4:6,2:4)  = -skew(r1);
    Aeq(4:6,5:7)  = -skew(r2);
    Aeq(4:6,8:10) = -skew(r3);
    beq(1:3) = [0;0;80];
    beq(4:6) = -skew(r1+r2+r3)*[0;0;20];
    
    x0 = ones(10, 1); % alpha, f1, f2, f3

    options = optimoptions('fmincon', 'SpecifyObjectiveGradient',true,...
      'SpecifyConstraintGradient',true,...
      'Display', 'none', 'Algorithm','interior-point') ;
    [xsol, ~, ~, ~, multiplier] = ...
      fmincon(@(x)low_cost(x), x0, [], [], Aeq, beq, [], [],...
      @(x)low_const(x), options);
    t = xsol(1);
    w(k) = exp(-beta * t);
    expsum = expsum + w(k);
    lambda(k, :) = multiplier.ineqnonlin;
    nu(k, :) = multiplier.eqlin;
    x(k, :) = xsol; % forces
  end
  obj = 1/beta * log(expsum); % maximization
  w = w / expsum;

  if nargout >1
    j = zeros(18, 1);
    for k = 1:N
      Q = zeros(19, 19);
      b = zeros(19, 9);

      yk = y(k, :);
      lambdak = lambda(k,:);
      xk = x(k, :);
      nu2 = nu(k, 4:6);
      f1 = xk(2:4);
      f2 = xk(5:7);
      f3 = xk(8:10);

      Q(1,17:19) = yk;
      Q(17:19,1) = yk';

      Q(2:4,2:4) = 2 * lambdak(1)*eye(3);
      Q(5:7,5:7) = 2 * lambdak(2)*eye(3);
      Q(8:10,8:10) = 2 * lambdak(3)*eye(3);

      Q(2:4,11)  = 2 * f1;
      Q(5:7,12)  = 2 * f2;
      Q(8:10,13) = 2 * f3;
      Q(11,2:4)  = 2 * lambdak(1)*f1;
      Q(12,5:7)  = 2 * lambdak(2)*f2;
      Q(13,8:10) = 2 * lambdak(3)*f3;

      Q(11, 11) = f1*f1'-900;
      Q(12, 12) = f2*f2'-900;
      Q(13, 13) = f3*f3'-900;

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

      b(2:4,1:3) = skew(nu2);
      b(5:7,4:6) = skew(nu2);
      b(8:10,7:9) = skew(nu2);
      b(17:19,1:3) = -skew(f1-[0,0,20]);
      b(17:19,4:6) = -skew(f2-[0,0,20]);
      b(17:19,7:9) = -skew(f3-[0,0,20]);

      dxdr = Q\b;
      dedr = dxdr(1,1:9);
      j(10:18) = j(10:18) + w(k) * dedr';
    end
    j = -j;
  end
end

function [obj, j] = low_cost(x)
  obj = -x(1);
  if nargout >1
    j = zeros(10, 1);
    j(1) = -1;
  end
end

function [g,h,gg,gh] = low_const(x)
  g(1) = x(2:4)'*x(2:4)-900;
  g(2) = x(5:7)'*x(5:7)-900;
  g(3) = x(8:10)'*x(8:10)-900;
  h = [];

  if nargout > 2
    gg = zeros(9,3);
    gg(2:4, 1) = 2*x(2:4);
    gg(5:7, 2) = 2*x(5:7);
    gg(8:10, 3) = 2*x(8:10);

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
