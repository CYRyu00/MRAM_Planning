% Newton-Euler Inverse Dynamics for n-DOF Robot Manipulator
function tau = newton_euler_ID(n, A, M, mass, inertia, F_tip, q, qd, qdd, g)
    % n: number of joints
    % A: A_i     R^6    screw axis for joint i, expressed in frame i , i = 1 ~ n
    % M: M_i,i-1 SE(3)  at home position, configuration of frame i-1 in frame i, i = 1 ~ n+1  
    % mass: mass of each link
    % inertia: inertia tensor of each link at that CoM
    % F_tip: wrench at the end-effector
    % q, qd, qdd: joint positions, velocities, accelerations
    % g: gravity vector

    %% Initialization
    T = cell(1, n+1); % T{i}  :  T_i,i-1 
    V = cell(1, n+1); % V{i}  : (w_i, v_i) spatial velocity, twist of link i, expressed in frame i 
    Vd= cell(1,n+1); % Vd{i} : (wd_i, vd_i)
    F = cell(1, n+1); % F{i}  : wrench transmitted through joint i to link frame i, expressed in frame i 
    G = cell(1, n);   % G{i}  : Spatial inertia matrix of link i, expressed in frame i 
    tau = zeros(n,1); % tau(i) : joint torque
    
    % Base conditions
    V_0 = [0;0;0;0;0;0];
    Vd_0 = [0;0;0;-g]; 
    F{n+1} = F_tip;
    T{n+1} = M{n+1};
    
    %% Forward Recursion (Kinematics)
    for i = 1:n
        T{i} = expm( -hat6(A{i}) *q(i)) *M{i};
        
        if i==1
            V{i}  = Ad_SE3(T{i}) *V_0 + A{i} *qd(i);
            Vd{i} = Ad_SE3(T{i}) *Vd_0 + ad_R6(V{i}) *A{i} *qd(i) + A{i} *qdd(i);
        else
            V{i}  = Ad_SE3(T{i}) *V{i-1} + A{i} *qd(i);
            Vd{i} = Ad_SE3(T{i}) *Vd{i-1} + ad_R6(V{i}) *A{i} *qd(i) + A{i} *qdd(i);
        end
        %disp(V{i})
        %disp(Vd{i})
    end

    %% Backward Recursion (Dynamics)
    for i = n:-1:1
        G{i} = [ inertia{i} , zeros(3,3) ; zeros(3,3) , mass{i}* eye(3,3) ];
        F{i} = Ad_SE3(T{i+1})' *F{i+1} + G{i} *Vd{i} - ad_R6(V{i})' *G{i} *V{i};
        tau(i) = F{i}' *A{i};
    end
end
