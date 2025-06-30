function R = rpy2rot(roll, pitch, yaw)
% Convert roll-pitch-yaw angles to SO(3) rotation matrix (extrinsic XYZ order)
% Input: roll, pitch, yaw in radians
% Output: 3x3 rotation matrix

% Rotation about X-axis (roll)
Rx = [1, 0, 0;
      0, cos(roll), -sin(roll);
      0, sin(roll), cos(roll)];

% Rotation about Y-axis (pitch)
Ry = [cos(pitch), 0, sin(pitch);
      0, 1, 0;
      -sin(pitch), 0, cos(pitch)];

% Rotation about Z-axis (yaw)
Rz = [cos(yaw), -sin(yaw), 0;
      sin(yaw), cos(yaw), 0;
      0, 0, 1];

% Extrinsic XYZ: R = Rz * Ry * Rx
R = Rz * Ry * Rx;
end
