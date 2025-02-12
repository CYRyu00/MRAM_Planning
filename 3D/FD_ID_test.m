%% Door : 4R robot
n = 4;
L1 = 0.2;

M = cell(1,n+1);
A = cell(1,n);
w = { [0;0;1],[0;0;1],[1;0;0],[0;1;0] };
r = { [1;2;3],[4;5;6],[7;8;9],[10;11;12] };
for i=1:n
    A{i} = [ w{i}; -hat(w{i})*r{i}];
end

%M{i} is M_i,i-1
M{1} = inv( [ eye(3,3), [1;2;3] ; 0,0,0,1 ]); 
M{2} = inv( [ eye(3,3), [4;5;6] ; 0,0,0,1 ]);
M{3} = inv( [ [0,1,0; -1,0,0; 0,0,1] , [7;8;9] ; 0,0,0,1 ]);
M{4} = inv( [ eye(3,3), [10;11;12] ; 0,0,0,1 ]);
M{5} = inv( [ eye(3,3), [13;14;15] ; 0,0,0,1 ]);


mass = {1,2,3,4};
inertia = {eye(3)*1, eye(3)*2, eye(3)*3, eye(3)*4};
g = [0;0;-9.81];


q = [pi/6; pi/5;pi/4;pi/3];
qd = [0.1; 0.2;0.3;0.4];
qdd = [0.5; 0.6;0.7;0.8];

F_tip = [1;2;3;4;5;6];
J = get_body_jacobian(n,A,M,q);
%%

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