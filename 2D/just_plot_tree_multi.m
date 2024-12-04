
x=x_opt;
u=u_opt;
L=L_arr;
num_up = worst_num_up;
robot = generate_model_multi_AMs(L,num_up);

figure('Position',[100,100,1000,600])
x0=[0,-pi/3,0,0]';
show(robot,x(1,1:2)'+x0(1:2));

view(2)
ax = gca;
ax.Projection = 'orthographic';
hold on

framesPerSecond = 1/dt*slow_factor;
rate = rateControl(framesPerSecond);


global params
mu = params(7);
r = params(8);
d = params(9);
A = [r r -r -r;-r r r -r;mu -mu mu -mu; zeros(2,4);ones(1,4)];

ax.View =[45,30];
axis([-0.7 0.7 -0.7 0.7 -1.3 1])
%%
for i = 1:N
    
    show(robot,x(i,1:2)'+ x0(1:2),'PreservePlot',false);
    drawnow
    hold on

    colors = ['r', 'g', 'b'];
    for j=1:length(num_up)
        tform1 = getTransform(robot, x(i,1:2)'+ x0(1:2), sprintf('AM%d-1', j));
        pos1 = tform1(1:3, 4);
        R1= tform1(1:3, 1:3);
       
        %thrust 1
        F_b = A*[u(i,2*j-1),0,0,0]';
        f_w = R1*F_b(4:6);   
        p1 = pos1 + R1*[r;r;0];
        p2 = p1 + f_w*scale;
        plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],colors(mod(j-1, 3) + 1), 'LineWidth',1);
    
        %thrust 2
        F_b = A*[0,u(i,2*j),0,0]';
        f_w = R1*F_b(4:6);   
        p1 = pos1 + R1*[-r;r;0];
        p2 = p1 + f_w*scale;
        plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],colors(mod(j-1, 3) + 1), 'LineWidth',1);
        
        if false
            hold off
            ax = gca;
            ax.Projection = 'orthographic';
            ax.View =[0,0];
            axis([-0.5 max(x_opt(:,1))+1.0 -0.5 0.5])
        end
    end

    waitfor(rate);
end
hold off