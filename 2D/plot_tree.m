function plot_tree(x_opt,u_opt, N_st,N_u,dt,N,slow_factor)
N_st = 4;N_u=2;
x=x_opt;
u=u_opt;
robot = generate_model_1AM();
%%

%close all

figure
x0=[0,-pi/2,0,0]';
show(robot,x(1,1:2)'+x0(1:2));

view(2)
ax = gca;
ax.Projection = 'orthographic';
hold on

framesPerSecond = 1/dt*slow_factor;
rate = rateControl(framesPerSecond);

scale=1.0;

global params
mu = params(7);
r = params(8);
A = [r r -r -r;-r r r -r;mu -mu mu -mu; zeros(2,4);ones(1,4)];

ax.View =[0,0];
axis([-0.5 5 -0.5 0.5])
for i = 1:N
    show(robot,x(i,1:2)'+ x0(1:2),'PreservePlot',false);
    drawnow
    hold on

    %AM1
    tform1 = getTransform(robot, x(i,1:2)'+ x0(1:2), 'AM1');
    pos1 = tform1(1:3, 4);
    R1= tform1(1:3, 1:3);
   
    %thrust 1
    F_b = A*[u(i,1),0,0,0]';
    f_w = R1*F_b(4:6);   
    p1 = pos1 + R1*[r;r;0];
    p2 = p1 + f_w*scale;
    plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'r', 'LineWidth',1);

    %thrust 2
    F_b = A*[0,u(i,2),0,0]';
    f_w = R1*F_b(4:6);   
    p1 = pos1 + R1*[-r;r;0];
    p2 = p1 + f_w*scale;
    plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'g', 'LineWidth',1);

    waitfor(rate);
end
hold off