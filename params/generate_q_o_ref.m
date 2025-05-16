function [q_o_ref] = generate_q_o_ref(x_0, x_f, N, dt, t0 ,t1)
N = N - t1 / dt *2;
T = N * dt;
q_o_ref = zeros(4, N + 1);
for j=1 : N + 1
    t = (j - 1) * dt;
    if t < t0
        q_o_ref(:, j) = [0 ; pi/6 * (cos(pi / t0 * t) - 1) / 2; - pi/6 * (cos(pi / t0 * t) - 1) / 2; pi/4 * (cos(pi / t0 * t) - 1) / 2 ];
        %q_o_ref(:, j) = [0 ; pi/6 * (cos(pi / t0 * t) - 1) / 2; 0; pi/4 * (cos(pi / t0 * t) - 1) / 2 ];
    elseif t < T-t0
        q_o_ref(:, j) = [(x_f(1) - x_0(1)) / (T - t0) * (t - t0) ; -pi/6; pi/6; -pi/4];
        %q_o_ref(:, j) = [(x_f(1) - x_0(1)) / (T - t0) * (t - t0) ; -pi/6; 0; -pi/4];
    else
        q_o_ref(:, j) = [(x_f(1) - x_0(1)) / (T - t0) * (t - t0) ; pi/6 * (cos(pi / t0 * t - pi / t0 * T) - 1) / 2; ...
                          -pi/6 * (cos(pi / t0 * t - pi / t0 * T) - 1) / 2; pi/4 * (cos(pi / t0 * t - pi / t0 * T) - 1) / 2];
        %q_o_ref(:, j) = [(x_f(1) - x_0(1)) / (T - t0) * (t - t0) ; pi/6 * (cos(pi / t0 * t - pi / t0 * T) - 1) / 2; ...
        %                  0; pi/4 * (cos(pi / t0 * t - pi / t0 * T) - 1) / 2];
    end
end
q_o_ref = [repmat(x_0(1:4), 1, t1 / dt), q_o_ref, repmat(x_f(1:4), 1, t1 / dt)];
end

