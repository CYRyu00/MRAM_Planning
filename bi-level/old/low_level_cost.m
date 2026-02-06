function [obj, j] = low_level_cost(z)
  y = z(1:3);
  mu = z(4:6);
  obj1 = -80*mu(3);

  r1 = [ 0  ;  0.8; 0];
  r2 = [-0.8; -0.4; 0];
  r3 = [ 0.8; -0.4; 0];
  a = skew([0,0,20]) * (r1+r2+r3);
  obj2 = a'*y';

  obj31 = norm(-skew(r1)*y' + mu');
  obj32 = norm(-skew(r2)*y' + mu');
  obj33 = norm(-skew(r3)*y' + mu');
  obj3 = obj31+obj32+obj33;
  obj = obj1 + obj2 + 30 * obj3;
  if nargout >1
    j1 = [0;0;0;0;0;-80];
    j2 = [a;0;0;0];

    j31m = ((-skew(r1)*y'+mu')/obj31)';
    j32m = ((-skew(r2)*y'+mu')/obj32)';
    j33m = ((-skew(r3)*y'+mu')/obj33)';
    j31y = -j31m*skew(r1);
    j32y = -j32m*skew(r2);
    j33y = -j33m*skew(r3);
    j3m = 30 * (j31m+j32m+j33m)';
    j3y = 30 * (j31y+j32y+j33y)';
    j3 = [j3y;j3m];
    j = j1+j2+j3;
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