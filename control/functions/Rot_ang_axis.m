function R = Rot_ang_axis(theta, u)
if abs(theta) < 1e-8
    R = eye(3);
else
    R = eye(3) + sin(theta) * S(u) + (1 - cos(theta)) * S2(u);
end
end
