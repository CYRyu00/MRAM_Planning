function S = hat(v)
    % HAT Map: Converts a 3D vector into a skew-symmetric matrix
    % Input:
    %   v - A 3x1 vector
    % Output:
    %   S - A 3x3 skew-symmetric matrix

    % Validate input
    if ~isvector(v) || length(v) ~= 3
        error('Input must be a 3-dimensional vector.');
    end

    % Extract vector components
    v1 = v(1);
    v2 = v(2);
    v3 = v(3);

    % Construct skew-symmetric matrix
    S = [  0   -v3   v2;
          v3     0  -v1;
         -v2   v1    0 ];
end
