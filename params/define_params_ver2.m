function params = define_params_ver2()
I0 = [0.001202, 0.0, 0.0002;
      0.0, 0.01056, 0.0;
      0.0002, 0.0, 0.0195];
thrust_limit = 20.0;
m0 = 1.8;

mu = 0.02;
r = 0.25; 
d = 0.4; % com to com
kt = 0.1; % Nm/rad at the handle

c_1 = 0.3; % Nms/rad
c_2 = 0.01;

mass_door = [10,1,0.05];
handle_factor = 0.95 * 9.81 * 0.0475;

inertia = {eye(3) * 1, eye(3) * 0.1, eye(3) * 0.1, zeros(3, 3), zeros(3, 3)};
r_i_ci = {[0.5; -0.02; 0.05], [-0.05; 0; 0.08], [0; 0; -0.05], [0; 0.05; 0], zeros(3, 1)};

n = 4;
dh = [0, 0, 0.95, 0;   % [alpha, a, d, theta]
      -pi/2, 0.9, 0, 0;
      0, -0.1, 0.23, pi/2;
      pi/2, 0, 0, -pi/2;
      pi/2, 0, 0, 0];
gravity = [0; 0; -9.81];

params = {m0, I0, mu, r, d, thrust_limit, kt, c_1, c_2, mass_door, handle_factor, inertia, r_i_ci, n, dh, gravity};
end