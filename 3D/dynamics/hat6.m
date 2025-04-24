function xi_hat = hat6(xi)
    % HAT6: Converts a 6D twist vector into a 4x4 matrix in se(3)
    % Input:
    %   xi - A 6x1 vector [omega; v] where omega and v are 3x1 vectors
    % Output:
    %   xi_hat - A 4x4 matrix representing the twist in se(3)
    
    % Validate input
    if ~isvector(xi) || length(xi) ~= 6
        error('Input must be a 6-dimensional vector.');
    end
    
    % Extract angular and linear components
    omega = xi(1:3);
    v = xi(4:6);
    
    % Compute the skew-symmetric matrix for omega
    omega_hat = hat(omega);
    
    % Construct the se(3) matrix
    xi_hat = [omega_hat, v;
              0, 0, 0, 0];
end