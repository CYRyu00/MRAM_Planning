function robot = generate_door_duo(n, dh, r_i_ci, d, g, shape, mass_list, inertia_list, do_view, q)
robot = rigidBodyTree('DataFormat','column','MaxNumBodies',n+1);
robot.Gravity = g';

% Example dimensions: [length, width, height]
box_size = { [1, 0.02, 2], [0.15, 0.05, 0.03],[0.04, 0.04, 0.1],[d*1.8, 0.04, d*0.9] ,[0.01,0.01,0.01]} ; 

% Build each link and joint
for i = 1:4
    % Create rigid body
    body = rigidBody(['link', num2str(i)]);
    
    % Create joint and define axis
    joint = rigidBodyJoint(['joint', num2str(i)], 'revolute');
    joint.JointAxis = [0;0;1];
    
    % Set fixed transform based on r_positions
    c1 = cos(dh(i,1));c2 = cos(dh(i,4)); 
    s1 = sin(dh(i,1));s2 = sin(dh(i,4));
    tform = [c2    -s2     0 dh(i,2); 
         c1*s2 c1*c2 -s1 -s1*dh(i,3); 
         s1*s2 c2*s1  c1 c1*dh(i,3);
         0 0 0 1]; 
    setFixedTransform(joint, tform);
    
    % Attach joint to body
    body.Joint = joint;
    
    % Define mass and inertia
    body.Mass = mass_list{i};
    body.CenterOfMass = r_i_ci{i}';
    body.Inertia = [inertia_list{i}(1,1), inertia_list{i}(2,2), inertia_list{i}(3,3), 0, 0, 0]; % Ixx, Iyy, Izz, Ixy, Iyz, Ixz
    
    % Add visual box to the body
    if i == 4
        [core_row, core_col] = find(shape == 2);
        [AMs_rows, AMs_cols] = find(shape ~= 0);
        for j =1:length(AMs_cols)
            r = [ (core_col - AMs_cols(j)) * -2 *d ; (core_row - AMs_rows(j)) *d ;0];% p_j,core
            r = -[r(1); -r_i_ci{i}(2); r(2)];
            
            addVisual(body, 'Box', box_size{i}, trvec2tform( r'))
        end
    else
        addVisual(body, 'Box', box_size{i}, trvec2tform(r_i_ci{i}'));
    end
    addVisual(body, 'Sphere', 0.03, trvec2tform(r_i_ci{i}'));
    
    % Attach body to robot
    if i == 1
        addBody(robot, body, robot.BaseName);
    else
        addBody(robot, body, ['link', num2str(i-1)]);
    end
end

fprintf("Robot is Generated Successfully\n")

if do_view ==1
% Visualize the robot 
close all
figure('Color', 'w');
ax = show(robot,q);
title('Generated robot model');

% Adjust view angle, zoom, and lighting
view(135, 20); % Set custom view angle
camzoom(0.7);  % Zoom in
axis equal;    % Maintain aspect ratio
light('Position', [1 1 1], 'Style', 'infinite'); % Add lighting for better visibility
grid on;      % Add grid for reference
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
end

end



