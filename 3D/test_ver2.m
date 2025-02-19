% Example test case for 2R Robot
n = 2;
L1=0.2; L2=0.3;
mass = { 1, 2};
inertia = {eye(3)*0.0, eye(3)*0.0};
r_i_ci ={ [L1;0;0] , [L2;0;0] };

g = [0; -9.81; 0];
q = [pi/4; pi/6];
qd = [1.0; 2.0];
qdd = [3.0; 4.0];
DH = [0, 0, 0, 0; 0, L1, 0, 0;  0, L2, 0, 0];

F_ext = [1;2;3;4;5;6]; % Force in end-effector local frame
%F_ext = zeros(6,1);

disp('Computed Joint Torques:');
tic
tau_ver2 = newton_euler_inverse_dynamics(n, DH, mass, inertia, r_i_ci, g, q, qd, qdd, F_ext);
toc
disp(tau_ver2);

% 2R robot 
A = { [0; 0; 1; 0; L1; 0] , [0; 0; 1; 0; L2; 0] };
M = { [ eye(3,3), [-L1 ; 0; 0] ; 0,0,0,1 ], [ eye(3,3), [-L2 ; 0; 0] ; 0,0,0,1 ], [ eye(4,4)] };

g = [0;-9.81;0];
J = get_body_jacobian(n,A,M,q);

disp("Exp. Inverse Dynmaics:")
tic
tau_ID = newton_euler_ID(n, A, M, mass, inertia, F_ext, q, qd, qdd, g);
toc
disp(tau_ID)

m1 = mass{1}; m2 = mass{2};

M_q = [ (m2*L2^2 + (m1+m2)*L1^2 +2*m2*L1*L2*cos(q(2)) ) , m2*L2*(L2+L1*cos(q(2))) ; m2*L2*( L2 + L1*cos(q(2))) ,m2*L2^2];
h_q_qd = [-2*m2*L1*L2*sin(q(2))*qd(1)*qd(2) - m2*L1*L2*sin(q(2))*qd(2)^2; m2*L2*L1*sin(q(2))*qd(1)^2];
G_q = [(m1+m2)*L1*cos(q(1))*9.81 + m2*L2*cos(q(1)+q(2))*9.81 ;m2*L2* cos(q(1)+q(2))*9.81 ];

disp("Analytic Inverse Dynmaics:")
tic
tau_analytic =  M_q *qdd + h_q_qd + G_q + J'*F_ext;
toc
disp(tau_analytic)

%%
tau = [1;2];

disp("Forward Dynamics")

tic
qdd_ver2 = FD_ver2(n, DH, mass, inertia, r_i_ci, g, q, qd, tau, F_ext)
toc

tic
qdd_FD = FD(n, A, M, mass, inertia, F_ext, q, qd, tau, g)
toc

M_q = [ (m2*L2^2 + (m1+m2)*L1^2 +2*m2*L1*L2*cos(q(2)) ) , m2*L2*(L2+L1*cos(q(2))) ; m2*L2*( L2 + L1*cos(q(2))) ,m2*L2^2];
h_q_qd = [-2*m2*L1*L2*sin(q(2))*qd(1)*qd(2) - m2*L1*L2*sin(q(2))*qd(2)^2; m2*L2*L1*sin(q(2))*qd(1)^2];
G_q = [(m1+m2)*L1*cos(q(1))*9.81 + m2*L2*cos(q(1)+q(2))*9.81 ;m2*L2* cos(q(1)+q(2))*9.81 ];

tic
qdd_Analytic = M_q\(tau - h_q_qd - G_q - J'*F_ext)
toc

