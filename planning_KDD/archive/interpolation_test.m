% --- 1. 변수 정의 ---
% 경로점 (위치): 시작 0, 끝 0.3. (1행 2열)
dt = 0.1; N = 10/dt;

waypoints = [0, 0.3]; 
timePoints = [0, N*dt]; 
t_samples = (0:N)*dt;
velocityBoundaryCondition = [0, 0]; 
[q, qd, qdd, pp] = cubicpolytraj(waypoints, timePoints, t_samples, ...
                                 'VelocityBoundaryCondition', velocityBoundaryCondition);
x_n_init = [q; qd];

figure;
subplot(3,1,1);
plot(t_samples, q);
title('Position');
ylabel('x (m)');

subplot(3,1,2);
plot(t_samples, qd);
title('Velocity');
ylabel('xdot (m/s)');

subplot(3,1,3);
plot(t_samples, qdd);
title('Acceleration');
xlabel('Time (s)');
ylabel('xddot (m/s^2)');
%% friction test
clear;

v = 0:0.01:1;
tanh_v = tanh(1e1*v.^2);

figure
plot(v, tanh_v)

