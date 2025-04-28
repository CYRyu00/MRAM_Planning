function tau = get_tau(x,u,l)
     
    g=9.8;
    global params
    m1 = params(1);
    m2 = params(2);
    lp = params(3);
    lg = params(4);
    m0 = params(5);
    I0 = params(6);
    mu = params(7);
    r = params(8);

    
    %Calculate Generalized force
    
    A = [0 0;-2*r 2*r; 0 0; 0 0; 0 0 ; 2 2];
    F_b = A*u;% [moment ; force]
    R = [-sin(x(2)) 0 cos(x(2)); 0 -1 0; cos(x(2)) 0 sin(x(2))];
    f_w = R*F_b(4:6);
    J = [0 0;0 1;0 0; 1 l*cos(x(2)); 0 0; 0 l*sin(x(2))];
    tau = J'*[F_b(1:3);f_w];
end