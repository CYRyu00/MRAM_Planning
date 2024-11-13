function  robot = generate_model_1AM()
    robot = rigidBodyTree('DataFormat','column','MaxNumBodies',5);
    
    global params
    m1 = params(1);
    m2 = params(2);
    lp = params(3);
    lg = params(4);
    m0 = params(5);
    I0 = params(6);
    mu = params(7);
    r = params(8);

    
    g=9.8;
    
    body1 = rigidBody('cart');
    body1.Mass = m1;
    body1.Inertia = zeros(1,6);
    joint1 = rigidBodyJoint('joint1', 'prismatic');
    setFixedTransform(joint1,trvec2tform([0 0 0]));
    joint1.JointAxis = [1 0 0];
    joint1.HomePosition = 0.2;
    body1.Joint = joint1;
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
    setFixedTransform(joint3, trvec2tform([lp, 0, 0]));
    body3.Joint = joint3;
    addBody(robot, body3, 'link');
    
    body4 = rigidBody('AM1');
    body4.Mass = m0;
    body4.Inertia = zeros(1,6);
    joint4 = rigidBodyJoint('fix2','fixed');
    tform4 = [diag([-1,-1,1]), [lg;0;0];0 0 0 1 ];
    setFixedTransform(joint4, tform4);
    body4.Joint = joint4;
    addBody(robot, body4, 'sphere');
    
    robot.Gravity = [0,0,-g];
    fprintf("\nCart-pole-AM model Generated\n")
