function points = fibonacci_sphere(N)
    % Generates N quasi-uniformly distributed points on a sphere using Fibonacci lattice.
    % Returns an Nx3 matrix where each row is a point (x, y, z).

    golden_angle = pi * (3 - sqrt(5)); % Golden angle in radians
    points = zeros(N, 3); % Preallocate for efficiency

    for i = 0:N-1
        x = 1 - 2 * i / (N - 1); % Linearly space y from 1 to -1
        radius = sqrt(1 - x^2);  % Compute the radius at height y
        theta = i * golden_angle; % Azimuthal angle

        y = radius * cos(theta);
        z = radius * sin(theta);

        points(i+1, :) = [x, y, z]; % Store the computed point
    end
end
