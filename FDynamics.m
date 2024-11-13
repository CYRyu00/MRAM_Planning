function xdot = FDynamics(x,u)
     
    global params
    m1 = params(1);
    m2 = params(2);
    lp = params(3);
    lg = params(4);
    m0 = params(5);
    I0 = params(6);
    mu = params(7);
    r = params(8);

    g=9.8;
    l1 = lp+lg;
    % q1 q2 q1dot q2dot
    M = [m1+m2+m0, (m2*lp+m0*l1)*cos(x(2));(m2*lp+m0*l1)*cos(x(2)) , m2*lp^2+m0*l1^2+I0];
    C = [0 , -(m2*lp+m0*l1)*sin(x(2))*x(4); 0 0];
    G = [0 ; (m2*lp+m0*l1)*g*sin(x(2))];
    
    %Calculate Generalized force
    A = [r r -r -r;-r r r -r;mu -mu mu -mu; zeros(2,4);ones(1,4)];
    F_b = A*u(1:4);% [moment ; force]
    R = [-sin(x(2)) 0 cos(x(2)); 0 -1 0; cos(x(2)) 0 sin(x(2))];
    f_w = R*F_b(4:6);
    J = [0 0;0 1;0 0; 1 l1*cos(x(2)); 0 0; 0 l1*sin(x(2))];
    tau = J'*[F_b(1:3);f_w];
    
    xdot = [x(3:4); M\(tau - C*x(3:4)-G)];

end