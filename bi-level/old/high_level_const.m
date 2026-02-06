function [g,h, gg, gh] = high_level_const(z)
  global sp;

  x1 = z(1:3, 1);
  x2 = z(4:6, 1);
  x3 = z(7:9, 1);
  r1 = z(10:12, 1);
  r2 = z(13:15, 1);
  r3 = z(16:18, 1);

  g(1) = x1(3) - 0.02;
  g(2) =-x1(3) - 0.02;
  g(3) = x2(3) - 0.02;
  g(4) =-x2(3) - 0.02;
  g(5) = x3(3) - 0.02;
  g(6) =-x3(3) - 0.02;

  h(1) = x1'*x1-1;
  h(2) = x2'*x2-1;
  h(3) = x3'*x3-1;
  h(4:6) = r1 - sp.surface(x1) - 0.3 * x1;
  h(7:9) = r2 - sp.surface(x2) - 0.3 * x2;
  h(10:12) = r3 - sp.surface(x3) - 0.3 * x3;
  if nargout >2
    gg = zeros(6, 18);
    gh = zeros(12, 18);

    gg(1, 3) = 1;
    gg(2, 3) =-1;
    gg(3, 6) = 1;
    gg(4, 6) =-1;
    gg(5, 9) = 1;
    gg(6, 9) =-1;
    gg = gg';

    gh(1, 1:3) = 2*x1;
    gh(2, 4:6) = 2*x2;
    gh(3, 7:9) = 2*x3;
    gh(4:6, 1:3) = -0.3 * eye(3) - sp.dsurface(x1);
    gh(4:6, 10:12) = eye(3);
    gh(7:9, 4:6) = -0.3 * eye(3) - sp.dsurface(x2);
    gh(7:9, 13:15) = eye(3);
    gh(10:12, 7:9) = -0.3 * eye(3) - sp.dsurface(x3);
    gh(10:12, 16:18) = eye(3);
    gh = gh';
  end
end

