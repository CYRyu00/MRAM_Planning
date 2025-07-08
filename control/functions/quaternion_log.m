function v = quaternion_log(q)
    % Input: q - 4x1 unit quaternion [w; x; y; z]
    % Output: v - 3x1 vector representing logarithm (rotation vector)

    w = q(1);
    vec = q(2:4);
    vec_norm = norm(vec);

    if vec_norm < 1e-10
        v = [0; 0; 0];
    else
        theta = acos(w);
        v = (theta / vec_norm) * vec * 2; % 2*log(q) = rotation vector
    end
end