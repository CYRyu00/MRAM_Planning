function robot = generate_door(n,M,w,r_positions, g, mass_list,inertia_list, do_view,q)
robot = rigidBodyTree('DataFormat','column','MaxNumBodies',n+1);
robot.Gravity = g';
R = cell(1,n);
P = cell(1,n);

for i=1:n
    M_  = M{i};
    R{i} = M_(1:3,1:3); % R_i,i-1
    P{i} = M_(1:3,4);   % P_i,i-1
end

% Example dimensions: [length, width, height]
box_size = { [1, 0.02, 2], [0.15, 0.05, 0.03],[0.08, 0.02, 0.015],[0.2, 0.4, 0.04] } ; 

% Build each link and joint
for i = 1:4


    % Create rigid body
    body = rigidBody(['link', num2str(i)]);
    
    % Create joint and define axis
    joint = rigidBodyJoint(['joint', num2str(i)], 'revolute');
    joint.JointAxis = w{i};
    
    % Set fixed transform based on r_positions
    if i==1
        tform = inv(M{i})*trvec2tform( r_positions{i}');
    else    
        tform = inv(trvec2tform( r_positions{i-1}'))*inv(M{i})*trvec2tform( r_positions{i}');
    end
    setFixedTransform(joint, tform);
    
    % Attach joint to body
    body.Joint = joint;
    
    % Define mass and inertia
    body.Mass = mass_list{i};
    body.CenterOfMass = r_positions{i}';
    body.Inertia = [inertia_list{i}(1,1), inertia_list{i}(2,2), inertia_list{i}(3,3), 0, 0, 0]; % Ixx, Iyy, Izz, Ixy, Iyz, Ixz
    
    % Add visual box to the body
    addVisual(body, 'Box', box_size{i},inv(trvec2tform( r_positions{i}')));
    addVisual(body, 'Sphere', 0.03, inv(trvec2tform( r_positions{i}')));

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
%q = [pi/3;-pi/6;-pi/6;-pi/12];
ax = show(robot,q);

% Adjust view angle, zoom, and lighting
title('Generated Robot Model');
view(135, 20); % Set custom view angle
camzoom(0.7);  % Zoom in
axis equal;    % Maintain aspect ratio
light('Position', [1 1 1], 'Style', 'infinite'); % Add lighting for better visibility
grid on;      % Add grid for reference
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
end