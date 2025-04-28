function x_interp = generate_x_interp(x_0, x_f, N, dt, t1)
    x_interp = zeros(N + 1, 8);
    vel_x1 = (x_f(1) - x_0(1)) / ( (N - 2 * t1) * dt); 
    vel_x2 = (x_f(2) - x_0(2)) / ( (N - 2 * t1) * dt);
    vel_x3 = (x_f(3) - x_0(3)) / ( (N - 2 * t1) * dt); 
    vel_x4 = (x_f(4) - x_0(4)) / ( (N - 2 * t1) * dt);
    for k = 1:8
        x_interp(t1 / dt : (end - 1 - t1 / dt), k) = linspace(x_0(k), x_f(k), N + 1 - 2 * t1 / dt)';
    end
    for k = t1 / dt : (N - t1 / dt)
        x_interp(k, 5:8) = [vel_x1; vel_x2; vel_x3; vel_x4];
end
end

