function  robot = generate_model_2AM(shape_)
    robot = rigidBodyTree('DataFormat','column','MaxNumBodies',5);
    m1=5; m2=2; mt=1.5; mq=1;
    l2=0.7; lg=0.3; d = 0.5;
    g=9.8;
    
    body1 = rigidBody('cart');
    body1.Mass = m1;
    body1.Inertia = zeros(1,6);
    joint1 = rigidBodyJoint('joint1', 'prismatic');
    setFixedTransform(joint1,trvec2tform([0 0 0]));
    joint1.JointAxis = [1 0 0];
    joint1.HomePosition = 0.2;
    body1.Joint = joint1;
    addVisual(body1,"Box",[0.2 0.2 0.2])
    addBody(robot, body1, 'base');
        
    body2 = rigidBody('link');
    body2.Mass = 0;
    body2.Inertia = zeros(1,6);
    joint2 = rigidBodyJoint('joint2','revolute');
    setFixedTransform(joint2, trvec2tform([0,0,0]));
    joint2.JointAxis = [0 -1 0];
    joint2.HomePosition = pi/2;
    body2.Joint = joint2;
    addBody(robot, body2, 'cart');
    
    body3 = rigidBody('sphere');
    body3.Mass = m2;
    body3.Inertia = zeros(1,6);
    joint3 = rigidBodyJoint('fix1','fixed');
    setFixedTransform(joint3, trvec2tform([l2, 0, 0]));
    body3.Joint = joint3;
    addVisual(body3,"Sphere",0.05)
    addBody(robot, body3, 'link');
    
    body4 = rigidBody('AM1');
    body4.Mass = mt;
    body4.Inertia = zeros(1,6);
    joint4 = rigidBodyJoint('fix2','fixed');
    tform4 = [diag([-1,-1,1]), [lg;0;0];0 0 0 1 ];
    setFixedTransform(joint4, tform4);
    body4.Joint = joint4;
    addVisual(body4,"Box",[0.3 0.3 0.1])
    addBody(robot, body4, 'sphere');
    
switch shape_
    case 1
        tform5 = trvec2tform([0, d, 0]);
    case 2
        tform5 = trvec2tform([-d, 0, 0]);
end

    body5 = rigidBody('AM2');
    body5.Mass = mq;
    body5.Inertia = zeros(1,6);
    addVisual(body5,"Box",[0.3 0.3 0.1])
    joint5 = rigidBodyJoint('fix3','fixed');
    
    setFixedTransform(joint5, tform5);
    body5.Joint = joint5;
    addBody(robot, body5, 'AM1');
    robot.Gravity = [0,0,-g];
    fprintf("\nCart-pole-AM model Generated\n")
