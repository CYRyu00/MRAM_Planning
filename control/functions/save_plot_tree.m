function save_plot_tree(robot,dh, params, x_opt,u_opt, dt,N,slow_factor, force_scale, shape)
x=x_opt;
u=u_opt;

n=4;
figure('Position',[50,100,800,600])

show(robot,x(1,1:4)');
view(2)
ax = gca;
hold on
dN = 10;
framesPerSecond = 1/dt*slow_factor/dN;
rate = rateControl(framesPerSecond);

video_filename = 'images/test.avi';
v = VideoWriter(video_filename); % Create a video writer object
v.FrameRate = framesPerSecond; % Set frame rate (adjust as needed)
open(v);

m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5}; thrust_limit= params{6};

ax.View =[150,25];
axis([-0.5 1.5 -0.5 1.5 0 2])
camzoom(1.2);


for i = 1:dN:N
    
    show(robot,x(i,1:4)','PreservePlot',false);
    drawnow
    hold on

    colors = ['r', 'g', 'b'];
    if true
        for j=1:1
        %T_04
        T_03 = getTransform(robot, x(i,1:4)', sprintf('link%d', 3));
        T_04 = getTransform(robot, x(i,1:4)', sprintf('link%d', 4));
        c1 = cos(dh(n+1,1)); c2 = cos(dh(n+1,4)); 
        s1 = sin(dh(n+1,1)); s2 = sin(dh(n+1,4));
        T_45= [  c2    -s2     0 dh(n+1,2); 
                 c1*s2 c1*c2 -s1 -s1*dh(n+1,3); 
                 s1*s2 c2*s1  c1 c1*dh(n+1,3);
                 0 0 0 1]; 
        T_05 = T_04*T_45;
        p_05 = T_05(1:3, 4);
        R_05 = T_05(1:3, 1:3);
        
        wrench = map_u2wrench_double( u(i,:)', shape , mu , r , d);
        
        % Force
        f_w = R_05 *wrench(4:6);   
        p1 = p_05;
        %p1 = p_05;
        p2 = p1 + f_w*force_scale;
        plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],colors( 1), 'LineWidth',1);
    
        % moment
        m_w = R_05 *wrench(1:3);   
        p1 = p_05;
        %p1 = p_05;
        p2 = p1 + m_w*force_scale;
        plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],colors(2), 'LineWidth',1);


        if false
            hold off
            ax = gca;
            ax.Projection = 'orthographic';
            ax.View =[0,0];
            axis([-0.5 max(x_opt(:,1))+1.0 -0.5 0.5])
        end
    end
    end
    frame = getframe(gcf);
    writeVideo(v, frame);
    waitfor(rate);
end
close(v);
disp(['Video saved to ', video_filename]);
hold off
end