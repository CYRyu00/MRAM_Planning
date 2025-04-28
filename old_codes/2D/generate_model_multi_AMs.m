function  robot = generate_model_multi_AMs(L,num_up)
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
    d = params(9);
    
    g=9.8;
    
    body1 = rigidBody('cart');
    body1.Mass = m1;
    body1.Inertia = zeros(1,6);
    joint1 = rigidBodyJoint('joint1', 'prismatic');
    setFixedTransform(joint1,trvec2tform([0 0 0]));
    joint1.JointAxis = [1 0 0];
    joint1.HomePosition = 0.2;
    body1.Joint = joint1;
    addVisual(body1,"Box",[0.1 0.1 0.1])
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
    addVisual(body3,"Sphere",0.03)
    addBody(robot, body3, 'link');
    

    for i=1:length(num_up)
        for j=1:num_up(i)
            body_name = sprintf('AM%d-%d', i,j);
            body4 = rigidBody(body_name);
            body4.Mass = m0;
            body4.Inertia = zeros(1,6);
            joint_name = sprintf('AM_fix%d-%d', i,j);
            joint4 = rigidBodyJoint(joint_name,'fixed');

            if j==1
            tform4 = [diag([-1,-1,1]), [L(i)-lp ;0;0];0 0 0 1 ];
            setFixedTransform(joint4, tform4);
            body4.Joint = joint4;
            addVisual(body4,"Box",[2*r 2*r r/2])
            addBody(robot, body4, 'sphere');
            else
            tform4 = [diag([1,1,1]), [0;-(j-1)*d;0];0 0 0 1 ];
            setFixedTransform(joint4, tform4);
            body4.Joint = joint4;
            addVisual(body4,"Box",[2*r 2*r r/2])
            addBody(robot, body4, sprintf('AM%d-1', i));
            end
        end
    end

    robot.Gravity = [0,0,-g];
    fprintf("\nCart-pole-AM model Generated\n")
