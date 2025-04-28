function J = get_body_jacobian(n, A, M, q)
    % n: number of joints
    % A: A_i     R^6   at home position, screw axis for joint i, expressed in frame i , i = 1 ~ n
    % M: M_i,i-1 SE(3) at home position, configuration of frame i-1 in frame i, i = 1 ~ n+1  
    B = cell(1, n);   % B{i}  :  at home position, screw axis for joint i, expressed in frame n+1 , i = 1 ~ n
    M_bi = eye(4,4); % M_b,b

    J = zeros(6,n);

    for i= n:-1:1
        M_bi = M_bi*M{i+1};
        B{i} = Ad_SE3(M_bi)*A{i};
        
        if i==n
            exp_B = eye(4,4);% exp(-B_n*q_n)*exp(-B_n-1*q_n-1) ...
        else
            exp_B = exp_B *expm(-hat6(B{i+1})*q(i+1));
        end

        J(:,i) = Ad_SE3(exp_B)*B{i};
        
    end

end
