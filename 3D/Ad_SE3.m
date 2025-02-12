function Ad_T = Ad_SE3(T)
    % ADJOINT_SE3: Computes the adjoint representation of a transformation matrix in SE(3)
    % Input:
    %   T - A 4x4 transformation matrix in SE(3)
    % Output:
    %   Ad_T - A 6x6 adjoint matrix

    % Validate input
    if ~isequal(size(T), [4, 4]) || abs(det(T(1:3,1:3)) - 1) > 1e-6
        error('Input must be a valid 4x4 transformation matrix in SE(3).');
    end

    % Extract rotation and translation
    R = T(1:3, 1:3);    % Rotation matrix
    p = T(1:3, 4);      % Translation vector

    % Compute skew-symmetric matrix of translation vector
    p_hat = hat(p);

    % Construct adjoint matrix
    Ad_T = [R, zeros(3,3);
            p_hat * R, R];
end