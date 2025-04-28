function tau_desired = get_tau_desired(x,xdot, L)
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


    term1 = m1 + m2;
    term2 = m2*lp;
    term3 = m2*lp^2;
    for i=1:length(L)
        term1 = term1 + m0;
        term2 = term2 + m0*L(i);
        term3 = term3 + m0*L(i)^2 + I0;
    end
    
   
    % q1 q2 q1dot q2dot
    M = [term1, term2*cos(x(2));term2*cos(x(2)) , term3];
    C = [0 , -term2*sin(x(2))*x(4); 0 0];
    G = [0 ; term2*g*sin(x(2))];
    
    tau_desired = M*xdot(3:4) + C*x(3:4) + G;
end