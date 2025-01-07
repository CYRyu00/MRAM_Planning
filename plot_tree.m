N_st = 4;N_u=4;
N = gridN;
if ~isempty(x_opt)
    x=x_opt;
    u=u_opt;

else
x = optimal(1:N*N_st);
u = optimal(N*N_st+1:end);
x = reshape(x,N,N_st);%row vectors
u = reshape(u,N,N_u);%row vectors
t = linspace(0,1,N)'*t_f;
dt = t_f/N;
end

robot = generate_model_1AM();
%%
x0=[0,-pi/2,0,0]';
figure(1)
show(robot,x(1,1:2)'+x0(1:2));
view(2)
ax = gca;
ax.Projection = 'orthographic';
ax.View =[45,30];
hold on
axis([-0.5 x_f(1)+0.5 -1 1])

framesPerSecond = 4;
rate = rateControl(framesPerSecond);

scale=0.1;

global params
mu = params(7);
r = params(8);
A1 = [r r -r -r;-r r r -r;mu -mu mu -mu; zeros(2,4);ones(1,4)];

for i = 1:gridN
    show(robot,x(i,1:2)'+ x0(1:2),'PreservePlot',false);
    drawnow
    hold on

    ax.View =[0,0];%[45,30];
    axis([-0.5 x_f(1)+0.5 -1 1])
    
    %AM1
    tform1 = getTransform(robot, x(i,1:2)'+ x0(1:2), 'AM1');
    pos1 = tform1(1:3, 4);
    R1= tform1(1:3, 1:3);
   
    %thrust 1
    F_b = A1*[u(i,1),0,0,0]';
    f_w = R1*F_b(4:6);   
    p1 = pos1 + R1*[r;r;0];
    p2 = p1 + f_w*scale;
    plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'r', 'LineWidth',1);

    %thrust 2
    F_b = A1*[0,u(i,2),0,0]';
    f_w = R1*F_b(4:6);   
    p1 = pos1 + R1*[-r;r;0];
    p2 = p1 + f_w*scale;
    plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'g', 'LineWidth',1);

    waitfor(rate);
end
%%
x0 = [0, -pi/2, 0, 0]';
figure(1)
show(robot, x(1, 1:2)' + x0(1:2));
view(2)
ax = gca;
ax.Projection = 'orthographic';
ax.View = [45, 30];
hold on
axis([-0.5 x_f(1) + 0.5 -1 1])

framesPerSecond = 1/dt;
rate = rateControl(framesPerSecond);

scale = 0.1;

global params
mu = params(7);
r = params(8);
A1 = [r r -r -r; -r r r -r; mu -mu mu -mu; zeros(2, 4); ones(1, 4)];

% Video Writer Setup
videoFileName = 'test.avi';
v = VideoWriter(videoFileName, 'Motion JPEG AVI');
v.FrameRate = framesPerSecond;
open(v);

for i = 1:gridN
    hold off
    show(robot, x(i, 1:2)' + x0(1:2), 'PreservePlot', false);
    drawnow
    hold on

    ax.View = [0, 0]; % [45,30];
    axis([-0.5 x_f(1) + 0.5 -1 1])
    
    % AM1
    tform1 = getTransform(robot, x(i, 1:2)' + x0(1:2), 'AM1');
    pos1 = tform1(1:3, 4);
    R1 = tform1(1:3, 1:3);
   
    % Thrust 1
    F_b = A1 * [u(i, 1), 0, 0, 0]';
    f_w = R1 * F_b(4:6);   
    p1 = pos1 + R1 * [r; r; 0];
    p2 = p1 + f_w * scale;
    plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'r', 'LineWidth', 1);

    % Thrust 2
    F_b = A1 * [0, u(i, 2), 0, 0]';
    f_w = R1 * F_b(4:6);   
    p1 = pos1 + R1 * [-r; r; 0];
    p2 = p1 + f_w * scale;
    plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'g', 'LineWidth', 1);

    % Capture the frame
    frame = getframe(gcf);
    writeVideo(v, frame);

    waitfor(rate);
end
hold off

% Close the video writer
close(v);
disp(['Video saved as ', videoFileName]);
%%

figure(1)
i= 39;
show(robot, x(i, 1:2)' + x0(1:2), 'PreservePlot', false);
ax = gca;
ax.Projection = 'orthographic';
hold on

ax.View = [0, 0]; % [45,30];
axis([-0.5 x_f(1) + 0.5 -1 1])


% AM1
tform1 = getTransform(robot, x(i, 1:2)' + x0(1:2), 'AM1');
pos1 = tform1(1:3, 4);
R1 = tform1(1:3, 1:3);

% Thrust 1
F_b = A1 * [u(i, 1), 0, 0, 0]';
f_w = R1 * F_b(4:6);   
p1 = pos1 + R1 * [r; r; 0];
p2 = p1 + f_w * scale;
plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'r', 'LineWidth', 1);

% Thrust 2
F_b = A1 * [0, u(i, 2), 0, 0]';
f_w = R1 * F_b(4:6);   
p1 = pos1 + R1 * [-r; r; 0];
p2 = p1 + f_w * scale;
plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'g', 'LineWidth', 1);
