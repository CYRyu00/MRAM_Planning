robot = rigidBodyTree("DataFormat","column");
base = robot.Base;
fixedbase = rigidBody("fixed_base");
arm1 = rigidBody("door");
arm2 = rigidBody("handle");
gripper = rigidBody("gripper");
AM = rigidBody("AM");

collBase = collisionCylinder(0.05,0.04); % cylinder: radius,length
collBase.Pose = trvec2tform([0 0 -0.04/2]);

coll1 = collisionBox(0.50,0.05,1.00); % box: length, width, height (x,y,z)
coll1.Pose = trvec2tform([0.05/2 0.50/2 1.00/2]);

coll2 = collisionBox(0.02,0.15,0.05); % box: length, width, height (x,y,z)
coll2.Pose = trvec2tform([0, -0.15/2, 0]);

%gripper
coll3 = collisionBox(0.10,0.05,0.02); % box: length, width, height (x,y,z)
coll3.Pose = trvec2tform([-0.10/2 0 0]);

dronebox = collisionBox(0.20,0.20,0.02);% sphere: radius
dronebox.Pose = trvec2tform([0 0 0]);

addCollision(fixedbase,collBase)
addCollision(arm1,coll1)
addCollision(arm2,coll2)
addCollision(gripper,coll3)
addCollision(AM,dronebox)

%Joint
jntBase = rigidBodyJoint("base_joint","fixed");
jnt1 = rigidBodyJoint("jnt1","revolute");
jnt2 = rigidBodyJoint("jnt2","revolute");
jnt3 = rigidBodyJoint("jnt3","revolute");
jntAM = rigidBodyJoint("AM_fixed_jnt","fixed");

jnt1.JointAxis = [0 0 1]; % z-axis
jnt2.JointAxis = [1 0 0];
jnt3.JointAxis = [1 0 0];
%jntGripper.JointAxis = [0 1 0]; % y-axis

jnt1.HomePosition = 0;
jnt2.HomePosition = 0;%-pi/6;
jnt3.HomePosition = 0;%pi/3;

setFixedTransform(jnt1,trvec2tform([0 0 0]))
setFixedTransform(jnt2,trvec2tform([0, 0 ,0]))
setFixedTransform(jnt3,trvec2tform([0,-0.15,0]))
setFixedTransform(jntAM,trvec2tform([(-0.10 -0.20/2) 0 0]))


bodies = {base,fixedbase,arm1,arm2,gripper,AM};
joints = {[],jntBase,jnt1,jnt2,jnt3,jntAM};

for i = 2:length(bodies) % Skip base. Iterate through adding bodies and joints.
            bodies{i}.Joint = joints{i};
            addBody(robot,bodies{i},bodies{i-1}.Name)
end
robot.Gravity =[0;0;-9.81];
%% show
x0=[0,-pi/2,0,0]';
show(robot,"Collisions","on");

view(2)
ax = gca;
ax.Projection = 'orthographic';
ax.View =[30,20];
