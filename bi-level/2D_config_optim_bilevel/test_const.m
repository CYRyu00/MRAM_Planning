function [g,h, gg, gh] = test_const(theta, num_AMs, f_min, f_max)
  % g<=0
  g = [f_min * ones(num_AMs, 1) - theta;
       theta - f_max * ones(num_AMs, 1)];
  % h=0
  h = [];
  if nargout > 2
    gg = zeros(num_AMs * 2, num_AMs);
    gh = [];

    gg(1:num_AMs, :) = -eye(num_AMs, num_AMs);
    gg(num_AMs+1:2*num_AMs, :) = eye(num_AMs, num_AMs);

    gg = gg';
  end
end

