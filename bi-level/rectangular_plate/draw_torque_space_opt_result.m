
u_ = [1, 0, 0];
or1 = zsol(10:12, 1); or2 = zsol(13:15, 1); or3 = zsol(16:18, 1);

Aeq = zeros(6,10);
beq = zeros(6,1);
Aeq(1:3,2:4) = eye(3);
Aeq(1:3,5:7) = eye(3);
Aeq(1:3,8:10) = eye(3);
Aeq(4:6,1) = u_';
Aeq(4:6,2:4) = -skew(or1);
Aeq(4:6,5:7) = -skew(or2);
Aeq(4:6,8:10) = -skew(or3);
beq(1:3) = 80*[0;0;1];
beq(4:6) = -20*(skew(or1+or2+or3)*[0;0;1]);
                                                       
theta1 = pi/180*(0:8:360);
theta2 = pi/180*(0:4:180);  
n1 = size(theta1,2);
X = zeros(n1,n1);
Y = zeros(n1,n1);
Z = zeros(n1,n1);
options = optimoptions('fmincon','SpecifyConstraintGradient',true,...
  'display', 'none', 'Algorithm','interior-point') ;
x0 = zeros(10, 1);
x0(4) = 20;
x0(7) = 20;
x0(10) = 20;
for i = 1:n1
  for j=1:n1
    x_ = cos(theta1(i))*sin(theta2(j));
    y_ = sin(theta1(i))*sin(theta2(j));
    z_ = cos(theta2(j));
    u_ = [x_, y_, z_];
    Aeq(4:6,1) = u_';
    [xsol, ~] = ...
      fmincon(@(x)low_level_cost(x), x0, [], [], Aeq, beq, [], [], @(x)low_level_const(x), options);
    alpha = -xsol(1);
    X(i, j) = alpha*x_;
    Y(i, j) = alpha*y_;
    Z(i, j) = alpha*z_;
  end
  i/n1
end

%%
s= mesh(X,Y,Z, 'FaceAlpha', 0.4, "EdgeAlpha",0.0);
s.FaceColor = 'flat';
s.FaceLighting = 'flat';
hold on
plot_sphere([0,0,0], radius);
plot3([0, -radius*dir_sol(1)], [0, -radius*dir_sol(2)], [0, -radius*dir_sol(3)], ...
  'k', "LineWidth", 2)
axis equal
xlabel("x")
ylabel("y")


function f = low_level_cost(x)
 f = -x(1);
end

function [g,h,gg,gh] = low_level_const(x)
  g(1) = x(2:4)'*x(2:4)-900;
  g(2) = x(5:7)'*x(5:7)-900;
  g(3) = x(8:10)'*x(8:10)-900;
  h = [];

  if nargout >2
    gg = zeros(10,3);
    gg(2:4,1) = 2*x(2:4);
    gg(5:7,2) = 2*x(5:7);
    gg(8:10,3) = 2*x(8:10);
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

function plot_sphere(c, r)
  theta1 = pi/180*(0:8:360);
  theta2 = pi/180*(0:4:180);  
  n1 = size(theta1,2);
  X = zeros(n1,n1);
  Y = zeros(n1,n1);
  Z = zeros(n1,n1);
  for i = 1:n1
    for j=1:n1
      x = c(1)+r*cos(theta1(i))*sin(theta2(j));
      y = c(2)+r*sin(theta1(i))*sin(theta2(j));
      z = c(3)+r*cos(theta2(j));
      X(i, j) = x;
      Y(i, j) = y;
      Z(i, j) = z;
    end
  end
  s= mesh(X,Y,Z, 'FaceAlpha', 0.5, "EdgeAlpha",0.5);
  s.FaceColor = 'flat';
  s.EdgeColor = 'interp';
  s.FaceLighting = 'flat';
end