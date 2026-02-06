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

