function At = get_At(x,L)
    
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
    A_f = [r r -r -r;-r r r -r;mu -mu mu -mu; zeros(2,4);ones(1,4)];
    R = [-sin(x(2)) 0 cos(x(2)); 0 -1 0; cos(x(2)) 0 sin(x(2))];
    At= [];
    for i= 1:length(L)
        J = [0 0;0 1;0 0; 1 L(i)*cos(x(2)); 0 0; 0 L(i)*sin(x(2))];
        A = J'*[eye(3,3),zeros(3,3);zeros(3,3),R]*A_f;
        At = [At , A];
    end
end