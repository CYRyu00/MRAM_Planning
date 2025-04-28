function xdot = FDynamics_multi(x,u)
     
    global params
    m1 = params(1);
    m2 = params(2);
    lp = params(3);
    lg = params(4);
    m0 = params(5);
    I0 = params(6);
    mu = params(7);
    r = params(8);
    global L
    num_AMs = length(L);
    term1 = m1 + m2;
    term2 = m2*lp;
    term3 = m2*lp^2;
    for i=1:num_AMs
        term1 = term1 + m0;
        term2 = term2 + m0*L(i);
        term3 = term3 + m0*L(i)^2 + I0;
    end
    % q1 q2 q1dot q2dot
    M = [term1, term2*cos(x(2));term2*cos(x(2)) , term3];
    C = [0 , -term2*sin(x(2))*x(4); 0 0];
    G = [0 ; term2*g*sin(x(2))];
    
    
    %Calculate Generalized force
    l1 = lp+lg;
    A = [r r -r -r;-r r r -r;mu -mu mu -mu; zeros(2,4);ones(1,4)];
    kron_A = kron(eye(num_AMs,num_AMs),A);
    F_b = kron_A*u(1:4);% [moment ; force]


    R = [-sin(x(2)) 0 cos(x(2)); 0 -1 0; cos(x(2)) 0 sin(x(2))];
    f_w = R*F_b(4:6);
    J = [0 0;0 1;0 0; 1 l1*cos(x(2)); 0 0; 0 l1*sin(x(2))];
    tau = J'*[F_b(1:3);f_w];
    
    xdot = [x(3:4); M\(tau - C*x(3:4)-G)];

end