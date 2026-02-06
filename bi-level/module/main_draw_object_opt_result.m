%% draw object
theta1 = pi/180*(0:4:360);
theta2 = pi/180*(0:2:180);  
n1 = size(theta1,2);
X = zeros(n1,n1);
Y = zeros(n1,n1);
Z = zeros(n1,n1);
for k = 1:object.size
  for i = 1:n1
    for j = 1:n1
      x_ = cos(theta1(i))*sin(theta2(j));
      y_ = sin(theta1(i))*sin(theta2(j));
      z_ = cos(theta2(j));
      s_ = object.surface(k, [x_;y_;z_]);
      X(i, j) = s_(1);
      Y(i, j) = s_(2);
      Z(i, j) = s_(3);
    end
  end
  s= mesh(X,Y,Z, 'FaceAlpha', 0.8, "EdgeAlpha",0.0);
  s.FaceColor = 'flat';
  s.FaceLighting = 'flat';
  hold on
end
grid on
axis equal
xlim([-1, 1])
ylim([-1, 1])
zlim([-1, 1])

%% draw drone & AM
s1 = object.surface(1, ox1); s2 = object.surface(2, ox2); s3 = object.surface(3, ox3);
ax1 = [s1(1), or1(1)]; ay1 = [s1(2), or1(2)]; az1 = [s1(3), or1(3)];
ax2 = [s2(1), or2(1)]; ay2 = [s2(2), or2(2)]; az2 = [s2(3), or2(3)];
ax3 = [s3(1), or3(1)]; ay3 = [s3(2), or3(2)]; az3 = [s3(3), or3(3)];

hold on
plot3(ax1, ay1, az1, "k", "linewidth", 2.0);
plot3(ax2, ay2, az3, "k", "linewidth", 2.0);
plot3(ax3, ay3, az2, "k", "linewidth", 2.0);