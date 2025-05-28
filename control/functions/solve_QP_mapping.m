function u = solve_QP_mapping(A, B, tau_des, f_des, u_min, u_max, R, eta)
% A: A_force
% B: A_tau
    eps_nz = 1.0e-10;
    
    % Moment matching (minimum-norm)
    B_pinv = pinv(B);
    u0     = B_pinv * tau_des;                      % particular solution
    Z      = null(B,'r');                         % basis of nullspace of B

    % If B has full row rank (rank 3) Z is 4×1;
    ksi = Z(:,1);

    % pre-compute quadratic coefficients  a₂ κ² + a₁ κ + a₀ = 0
    a2 = ksi' * A' * A * ksi;
    a1 = 2 * (u0'* A' * A * ksi);
    a0 = u0'* A' * A * u0 - f_des^2;
    fprintf("a2, a1, a0: %.2f, %.2f, %.2f\n", a2, a1, a0)

    roots_kappa = roots([a2 a1 a0]);              % at most two real roots
    fprintf("roots: %.2f, %.2f\n", roots_kappa(1), roots_kappa(2))
    u_feas  = [];                             % feasible thrust set
    for kappa = roots_kappa(:).'
        if abs(imag(kappa)) < eps_nz              % discard complex solutions
            u_cand = u0 + real(kappa) * ksi;
            u_feas = [u_feas u_cand];    
        end
    end

    % choose the feasible candidate with minimal Euclidean norm
    if ~isempty(u_feas)
        diff_min = inf;
        for idx = 1:length(u_feas(1, :))
            %disp(idx)
            f = A * u_feas(:, idx);
            diff = norm( f' * R * eta- f_des);
            %disp(A * u_feas(:, idx))
            if diff < diff_min
                u = u_feas(:, idx);
                diff_min = diff;
                min_idx = idx;
            end
        end
        %f = A * u_feas(:, min_idx);
        %disp(f' * R * eta- f_des)
        %disp(diff_min)
    else
        u = min(max(u0,u_min),u_max); 
    end
end
