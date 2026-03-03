function [obj, j] = high_cost_sparse(theta, beta, U, num_AMs, num_points, thrust_min, thrust_max, R_shape, r_cj, B_tau, B_lambda)
    expsum = 0;
    dim = 3;
    w = zeros(num_points, 1);
    % lambda = zeros(num_points, 3);
    nu = zeros(num_points, dim);
    x = zeros(num_points, 4*num_AMs+1);

    B_theta_3D = zeros(dim, 4, num_AMs);
    dB_dt_3D   = zeros(dim, 4, num_AMs);
    for i = 1:num_AMs
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
        'Algorithm', 'interior-point', ...
        'Display', 'none');
    lb = [-inf; thrust_min * ones(4*num_AMs, 1)];
    ub = [ inf; thrust_max * ones(4*num_AMs, 1)];

    h = zeros(4*num_AMs+1, 1);
    h(1) = -1;
    
    h_large = repmat(h, num_points, 1);
    beq_large = zeros(dim * num_points, 1);
    lb_large = repmat(lb, num_points, 1);
    ub_large = repmat(ub, num_points, 1);

    Aeq_cell = cell(num_points, 1);
    for k = 1:num_points
        Aeq_k = [U(:, k), -B_theta];
        Aeq_cell{k} = sparse(Aeq_k);
    end
    Aeq_large = sparse(blkdiag(Aeq_cell{:}));

    [xsol, ~, ~, ~, multiplier] = linprog(h_large, [], [], Aeq_large, beq_large, lb_large, ub_large, options);
    
    x_mat = reshape(xsol, 4*num_AMs+1, num_points);
    rho_all = x_mat(1, :);

    nu_mat = reshape(multiplier.eqlin, dim, num_points);
    w = exp(-beta * rho_all);
    expsum = sum(w);
    x = x_mat';
    nu = nu_mat';

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