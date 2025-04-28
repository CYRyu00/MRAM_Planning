function K_arr = compute_LQR_gains(A_arr, B_arr, Q, Qf, R, dt)
    N_sim = numel(A_arr);
    [nx, ~] = size(A_arr{1});
    [~, nu] = size(B_arr{1});
    P  = Qf;
    K_arr  = cell(1, N_sim);
    Ad_arr = cell(1, N_sim);
    Bd_arr = cell(1, N_sim);
    
    for i = 1:N_sim
        A = A_arr{i};
        B = B_arr{i};
        M = [A, B; zeros(nu, nx+nu)];
        Md = expm(M * dt);
        Ad_arr{i} = Md(1:nx,      1:nx);
        Bd_arr{i} = Md(1:nx, nx+1:end);
    end
    
    for i = N_sim:-1:1
        A = Ad_arr{i};
        B = Bd_arr{i};
        S = R + B' * P * B;
        K = S \ (B' * P * A);
        K_arr{i} = K;
        P = Q + A' * P * A - A' * P * B * (S \ (B' * P * A));
    end
end
