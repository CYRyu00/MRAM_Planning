%% GPT
n = 3;
DH_params = [0, 0.3, 0.2, pi/2;
             0, 0, 0.4, 0;
             0, 0, 0.1, -pi/2];
mass = [5, 3, 2];
com = {[0; 0; 0.15], [0.2; 0; 0], [0.1; 0; 0]};
inertia = {eye(3)*0.01, eye(3)*0.008, eye(3)*0.005};
q = [pi/4; pi/3; pi/6];
qd = [0.5; 0.4; 0.3];
qdd = [0.1; 0.2; 0.1];
g = [0; 0; 9.81];

tau = newton_euler_inverse_dynamics(n, DH_params, mass, com, inertia, q, qd, qdd, g);
disp(tau);
%% R-P example
n = 2;
L1 = 0.2;

A = { [0; 0; 1; 0; L1; 0] , [0; 0; 0; 1; 0; 0] };
M = { [ eye(3,3), [-L1 ; 0; 0] ; 0,0,0,1 ],  [ eye(3,3), [L1 ; 0; 0] ; 0,0,0,1 ], [ eye(4,4)] };

mass = {3, 4};
inertia = {eye(3)*0.02, eye(3)*0.01};
g = [0;-9.81;0];


q = [pi/6; 0.5];
qd = [0.0; 0.0];
qdd = [0.0; 0.0];

F_tip = [1;2;3;4;5;6];
J = get_body_jacobian(n,A,M,q);

tau_ID = newton_euler_ID(n, A, M, mass, inertia, F_tip, q, qd, qdd, g);

M_q = [ inertia{1}(3,3) + inertia{2}(3,3) + mass{1} *L1^2 + mass{2} *q(2)^2 , 0 ; 0 , mass{2} ];
C_q_qd = [ mass{2}*q(2)*qd(2) , mass{2}*q(2)*qd(1); -mass{2}*q(2)*qd(1),0];
G_q = [(mass{1}*L1 + mass{2}*q(2)) *9.81*cos(q(1)) ; mass{2} *9.81*sin(q(1)) ];
tau_analytic = M_q *qdd + C_q_qd*qd + G_q +J'*F_tip;

% Forward dynamics

q = [pi/6; 0.5];
qd = [0.0; 0.0];
F_tip = [1;2;3;4;5;6];

J = get_body_jacobian(n,A,M,q);
tau = [1;2];

qdd_FD = FD(n, A, M, mass, inertia, F_tip, q, qd, tau, g);

M_q = [ inertia{1}(3,3) + inertia{2}(3,3) + mass{1} *L1^2 + mass{2} *q(2)^2 , 0 ; 0 , mass{2} ];
C_q_qd = [ mass{2}*q(2)*qd(2) , mass{2}*q(2)*qd(1); -mass{2}*q(2)*qd(1),0];
G_q = [(mass{1}*L1 + mass{2}*q(2)) *9.81*cos(q(1)) ; mass{2} *9.81*sin(q(1)) ];
qdd_Analytic = M_q\(tau - C_q_qd*qd - G_q - J'*F_tip);

%% 2R robot 
n = 2;
L1 = 0.2; L2 = 0.3;

A = { [0; 0; 1; 0; L1; 0] , [0; 0; 1; 0; L2; 0] };
M = { [ eye(3,3), [-L1 ; 0; 0] ; 0,0,0,1 ], [ eye(3,3), [-L2 ; 0; 0] ; 0,0,0,1 ], [ eye(4,4)] };

mass = {3, 5};
inertia = {eye(3)*0.0, eye(3)*0.0};
g = [0;-9.81;0];

%inverse dynamics
q = [pi/6; pi/4];
qd = [0.2; 0.1];
qdd = [0.2; 0.3];

F_tip = [1;2;3;4;5;6];
J = get_body_jacobian(n,A,M,q);
disp("Inverse Dynmaics:")
tic
tau = newton_euler_ID(n, A, M, mass, inertia, F_tip, q, qd, qdd, g);
toc

tic
m1 = mass{1}; m2 = mass{2};
tau1= (m1+m2)*L1*cos(q(1))*9.81 + m2*L2*cos(q(1)+q(2))*9.81 ...
                + (m2*L2^2 + (m1+m2)*L1^2 +2*m2*L1*L2*cos(q(2)) )*qdd(1) ...
                + m2*L2*( L2 + L1*cos(q(2))) *qdd(2) ...
                -2*m2*L1*L2*sin(q(2))*qd(1)*qd(2) - m2*L1*L2*sin(q(2))*qd(2)^2;
tau2=  m2*L2*( cos(q(1)+q(2))*9.81 + L1*sin(q(2))*qd(1)^2+(L2+L1*cos(q(2)))*qdd(1) + L2*qdd(2) ) ;

M_q = [ (m2*L2^2 + (m1+m2)*L1^2 +2*m2*L1*L2*cos(q(2)) ) , m2*L2*(L2+L1*cos(q(2))) ; m2*L2*( L2 + L1*cos(q(2))) ,m2*L2^2];
h_q_qd = [-2*m2*L1*L2*sin(q(2))*qd(1)*qd(2) - m2*L1*L2*sin(q(2))*qd(2)^2; m2*L2*L1*sin(q(2))*qd(1)^2];
G_q = [(m1+m2)*L1*cos(q(1))*9.81 + m2*L2*cos(q(1)+q(2))*9.81 ;m2*L2* cos(q(1)+q(2))*9.81 ];

tau_analytic =  M_q *qdd + h_q_qd + G_q + J'*F_tip;
toc

%forward dynamicss

q = [pi/6; pi/4];
qd = [0.2; 0.1];

F_tip = [1;2;3;4;5;6];
J = get_body_jacobian(n,A,M,q);
tau = [1;2];

disp("Forward Dynamics")
tic
qdd_FD = FD(n, A, M, mass, inertia, F_tip, q, qd, tau, g);
toc
tic
M_q = [ (m2*L2^2 + (m1+m2)*L1^2 +2*m2*L1*L2*cos(q(2)) ) , m2*L2*(L2+L1*cos(q(2))) ; m2*L2*( L2 + L1*cos(q(2))) ,m2*L2^2];
h_q_qd = [-2*m2*L1*L2*sin(q(2))*qd(1)*qd(2) - m2*L1*L2*sin(q(2))*qd(2)^2; m2*L2*L1*sin(q(2))*qd(1)^2];
G_q = [(m1+m2)*L1*cos(q(1))*9.81 + m2*L2*cos(q(1)+q(2))*9.81 ;m2*L2* cos(q(1)+q(2))*9.81 ];

qdd_Analytic = M_q\(tau - h_q_qd - G_q - J'*F_tip);
toc