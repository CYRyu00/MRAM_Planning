function ad_V = ad_R6(V)
    % Lie bracket of V1 and V2: ad_V1 * V2
    % Input:
    %   V - A 6x1 vector [omega; v] where omega and v are 3x1 vectors
    % Output:
    %   ad_V - A 6x6 adjoint matrix

    % Validate input
    if ~isvector(V) || length(V) ~= 6
        error('Input must be a 6-dimensional vector.');
    end

     % Extract angular and linear components
    omega = V(1:3);
    v = V(4:6);
    
    % Compute the skew-symmetric matrix for omega
    omega_hat = hat(omega);
    v_hat = hat(v);
    % Construct the se(3) matrix
    ad_V = [omega_hat, zeros(3,3);
              v_hat , omega_hat];
end