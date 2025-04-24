function [x_dot_func, A_func, B_func] = define_dynamics(n, dh, shape, num_AMs, gravity, params)
    import casadi.*
    
    delta_inertia = MX.sym('delta_inertia',1,1);
    delta_k = MX.sym('delta_k',1,1);
    disturb = MX.sym('disturb',n,1);
    params = define_params();
    m0 = params{1} *delta_inertia; I0 = params{2}*delta_inertia; mass_door = params{10} *delta_inertia;
    mu = params{3}; r= params{4}; d= params{5}; thrust_limit= params{6}; kt=params{7} *delta_k; c_1=params{8}; c_2=params{9}; handle_factor = params{11};
    
    [AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);
    mass =  {mass_door(1), mass_door(2), mass_door(3), AM_mass};
    inertia = {eye(3)*1, eye(3)*0.1, eye(3)*0.1, AM_inertia, zeros(3,3)};
    r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};
    
    nx = n*2;
    nu = num_AMs*4;
    x = MX.sym('x', nx, 1); % [q; qd];
    u = MX.sym('u', nu, 1);
    tau = [-c_1*x(5);(-c_2*x(6) -kt*x(2) + mass{2}*handle_factor); 0; 0] + disturb;
    F_ext = map_u2wrench( u, shape , mu , r , d);
    
    qdd = FD_ver2(n, dh, mass, inertia, r_i_ci, gravity, x(1:4), x(5:8), tau, F_ext);
    
    x_dot = [x(5:8);qdd];
    
    A = jacobian(x_dot, x);
    B = jacobian(x_dot, u);
    
    x_dot_func = Function('x_dot_func', {x ,u ,delta_inertia, delta_k, disturb }, {x_dot});
    A_func = Function('A_func', {x, u, delta_inertia ,delta_k,disturb }, {A});
    B_func = Function('B_func', {x, u, delta_inertia, delta_k,disturb }, {B});
end