function [g,h,gg,gh] = high_level_constraint(u)
  g =[];
  h = u*u'-1;
  if nargout >2
    gg = [];
    gh = 2*u;
  end
end

