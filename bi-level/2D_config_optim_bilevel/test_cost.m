function [obj, j] = test_cost(theta, num_AMs, num_points)
  expsum = 0;
  beta = 5.0;
  w = zeros(num_points, 1);
  lambda = zeros(num_points, 3);
  nu = zeros(num_points, 6);
  x = zeros(num_points, 10);
  parfor k = 1:num_points    
    t = theta(1) * k; %TODO, primal setting
    w(k) = exp(-beta * t);
    expsum = expsum + w(k);
  end
  obj = 1/beta * log(expsum); % maximization
  w = w / expsum;

  if nargout > 1
    j = zeros(num_AMs, 1);
    for k = 1:num_points
        dr_dt = zeros(num_AMs, 1);
        dr_dt(1) = k;
        j = j + w(k) * dr_dt;
    end
    j = -j;
  end
end