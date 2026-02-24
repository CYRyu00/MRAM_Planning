function [obj, j] = high_cost(theta, U, num_AMs, num_points, thrust_min, thrust_max, R_shape, r_cj, B_tau, B_lambda)
    expsum = 0;
    beta = 1e-1;
    dim = 3;
    w = zeros(num_points, 1);
    % lambda = zeros(num_points, 3);
    nu = zeros(num_points, dim);
    x = zeros(num_points, 4*num_AMs+1);

    B_theta_3D = zeros(dim, 4, num_AMs);
    dB_dt_3D   = zeros(dim, 4, num_AMs);
    parfor i = 1:num_AMs
        c = cos(theta(i)); s = sin(theta(i));

        R = R_shape{i}' * [c, 0, -s; 0, 1, 0; s, 0, c]; % R_y(-theta)
        B_i = [R, hat(r_cj{i}) * R; zeros(3, 3), R] * [B_tau; zeros(2, 4); B_lambda];
        B_theta_3D(:, :, i) = [B_i(2, :); B_i(4, :); B_i(6, :)]; % TODO : 2D

        dR_dti = R_shape{i}' * [-s, 0, -c; 0, 0, 0; c, 0, -s];
        dB_dti = [dR_dti, hat(r_cj{i}) * dR_dti; zeros(3, 3), dR_dti] * [B_tau; zeros(2, 4); B_lambda];
        dB_dt_3D(:, :, i) = [dB_dti(2, :); dB_dti(4, :); dB_dti(6, :)];
    end
    B_theta = reshape(B_theta_3D, dim, 4*num_AMs);
    dB_dt = reshape(dB_dt_3D, dim, 4*num_AMs);
    
    options = optimoptions('linprog', ...
        'Algorithm', 'dual-simplex', ...
        'Display', 'none');
    lb = [-inf; thrust_min * ones(4*num_AMs, 1)];
    ub = [ inf; thrust_max * ones(4*num_AMs, 1)];
    parfor k = 1:num_points %TODO LP solver 개선 SPARSE, LP 전용 솔버 사용
        u_k = U(:, k);
        h = zeros(4*num_AMs+1, 1);
        h(1) = -1;
        Aeq = [u_k, -B_theta];
        beq = zeros(dim, 1);
        
        [xsol, ~, ~, ~, multiplier] = linprog(h, [], [], Aeq, beq, lb, ub, options);

        rho_k = xsol(1);
        % lambda(k, :) = multiplier.ineqnonlin;
        nu(k, :) = multiplier.eqlin;
        w(k) = exp(-beta * rho_k);
        x(k, :) = xsol;
        expsum = expsum + w(k);
    end
    obj = 1/beta * log(expsum); % maximization
    w = w / expsum; % dr_drho
    
    if nargout > 1
        j = zeros(num_AMs, 1);
        parfor k = 1:num_points
            dr_dt = zeros(num_AMs, 1); % drho_k / dtheta
            for i = 1:num_AMs
                e_i = zeros(num_AMs, 1);
                e_i(i) = 1;
                dr_dt = dr_dt + (nu(k, :) * dB_dt(:, 4*i-3:4*i) * x(k, 4*i-2:4*i+1)') * e_i;
            end
            j = j + w(k) * dr_dt;
        end
        j = -j;
    end
end