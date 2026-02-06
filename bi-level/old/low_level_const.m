function [g,h, gg, gh] = low_level_const(z)
  g = [];
  y = z(1:3);
  h = y*y'-1;
  if nargout >2
    gg = [];
    gh = [2*y, 0, 0, 0]';
  end
end

