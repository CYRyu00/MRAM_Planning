%% Door : 4R robot
n = 4;
% dynamic parameters of each module a
params = define_params();
m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5};

m=8;
shapes = generate_all_shapes(m);
shape = shapes{8}{47};
[AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);

M = cell(1,n+1);
A = cell(1,n);
w = { [0;0;1], [0;1;0], [1;0;0], [0;1;0] };
r = { [-0.5;-0.05;-1], [0.05;-0.07;0.0], [0.03;0;0], -AM_com };
for i=1:n
    A{i} = [ w{i}; -hat(w{i})*r{i}];
end

%M{i} is M_i,i-1
M{1} = inv( [ eye(3,3), [0.5;0.05;1] ; 0,0,0,1 ]); % base frame to door 
M{2} = inv( [ eye(3,3), [0.4;0.1;-0.05] ; 0,0,0,1 ]); % door to handle&gripper
M{3} = inv( [ [0,1,0; -1,0,0; 0,0,1] , [-0.03; 0.1 ; 0.0] ; 0,0,0,1 ]); % handle&gripper to link
M{4} = inv( [ eye(3,3), AM_com+[-0.05;0;0] ; 0,0,0,1 ]);% link to AM
M{5} = inv( [ eye(3,3), -AM_com ; 0,0,0,1 ]);% AM to core module

mass = {20, 3, 0.5, AM_mass};
inertia = {eye(3)*1, eye(3)*2, eye(3)*3, AM_inertia};
g = [0;0;-9.81];

do_show = 1; q = [0;0;0;0];
robot = generate_door(n,M,w,r, g, mass,inertia, do_show, q);
%% ID
q = [pi/6; pi/5;pi/4;-pi/3];
qd = [0.1; 0.2;0.3;0.4];
qdd = [0.5; 0.6;0.7;0.8];
F_tip = [1;2;3;4;5;6];
J = get_body_jacobian(n,A,M,q);
tau_ID = newton_euler_ID(n, A, M, mass, inertia, F_tip, q, qd, qdd, g)

% Forward dynamics
q = [pi/6; pi/5;pi/4;-pi/3];
qd = [0.1; 0.2;0.3;0.4];
F_tip = [1;2;3;4;5;6];

tau = [1;2;3;4];
J = get_body_jacobian(n,A,M,q);

qdd_FD = FD(n, A, M, mass, inertia, F_tip, q, qd, tau, g)