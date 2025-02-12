% Create a rigid body tree
robot = rigidBodyTree('DataFormat', 'column', 'MaxNumBodies', 5);

w = 0.7; h = 2.0;t=0.05;
% Door Panel
door = rigidBody('door');
door.Mass = 5; % Mass of the door
door.Inertia = [0.1 0.1 0.1 0 0 0]; % Simplified inertia
joint1 = rigidBodyJoint('hinge1', 'revolute');
setFixedTransform(joint1, trvec2tform([w/4, -t/2, h/2])); % Set the hinge at the base corner
joint1.JointAxis = [0 0 1]; % Rotate around Z-axis
door.Joint = joint1;
addVisual(door, "Box", [ w/2, t, h]); % Represent the door as a box
addBody(robot, door, 'base'); % Attach door to the base

% Door Handle
handle = rigidBody('handle');
handle.Mass = 0.5; % Mass of the handle
handle.Inertia = [0.01 0.01 0.01 0 0 0]; % Simplified inertia
joint2 = rigidBodyJoint('handleJoint', 'revolute');
setFixedTransform(joint2, trvec2tform([0.1*w, t, 0])); % Position the handle on the door
joint2.JointAxis = [1 0 0]; % Rotate around X-axis
handle.Joint = joint2;
addVisual(handle, "Box", [0.15, 0.05, 0.05]); % Represent the handle as a small box
addBody(robot, handle, 'door'); % Attach handle to the door

%% Visualize the Model
figure;
show(robot);
title('Door Model with Hinge and Handle');
