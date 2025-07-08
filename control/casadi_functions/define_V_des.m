function [V_des_func, grad_func, hess_func] = define_V_des()
    import casadi.*
    
    %V_des = MX.sym('V_des', 1, 1);
    q = MX.sym('q', 6 ,1);
    q_des = MX.sym('q', 6 ,1);
    m = MX.sym('m');
    g = MX.sym('g');
    K_X = MX.sym('K_X', 3, 1);
    K_r = MX.sym('K_r', 3, 3);

    v = q(1:3);
    theta = norm_2(v);
    if_else = @(cond, a, b) cond * a + (1-cond) * b;
    u = if_else(theta < 1e-10, [1;0;0], v/theta);
    
    ux = [  0    -u(3)  u(2);
           u(3)   0    -u(1);
          -u(2)  u(1)   0   ];
    
    R = if_else(theta < 1e-10, eye(3), ...
        eye(3) + sin(theta)*ux + (1-cos(theta))*(ux*ux));
   
    P_X = [q(4) * ( m * g + K_X(1) * (q_des(4) - q(4)) + K_X(2) * (q_des(5) - q(5)) + K_X(3) * (q_des(6) - q(6)) + 0.5 * q(4));
           q(5) * ( m * g + K_X(1) * (q_des(4) - q(4)) + K_X(2) * (q_des(5) - q(5)) + K_X(3) * (q_des(6) - q(6)) + 0.5 * q(5));
           q(6) * ( m * g + K_X(1) * (q_des(4) - q(4)) + K_X(2) * (q_des(5) - q(5)) + K_X(3) * (q_des(6) - q(6)) + 0.5 * q(6))];

    V_des = 0.5 * (q(1:3) - q_des(1:3))' * K_r * (q(1:3) - q_des(1:3)) + ...
            m* g * q(6) - [0, 0, 1] * R * P_X;

    grad_V = jacobian(V_des, q);
    hess_V = hessian(V_des, q);
    
    V_des_func = Function('V_des', {q, q_des, m, g, K_X, K_r}, {V_des});
    grad_func = Function('grad_func', {q, q_des, m, g, K_X, K_r}, {grad_V});
    hess_func = Function('hess_func', {q, q_des, m, g, K_X, K_r}, {hess_V});
end

