function [x_d_interp, u_d_interp] = interpolate_traj(x_d, u_d, t_plan, t_sim, do_plot)

    x_d_interp = interp1(t_plan, x_d, t_sim, 'pchip');  % 1001x8
    u_d_interp = interp1(t_plan(1:end-1), u_d, t_sim(1:end-1), 'pchip');  % 1000x40
    
    if do_plot
        figure('Position',[100 500 600 400])
        subplot(2,1,1)
        plot(t_plan, x_d, '.');
        title( {'$x_d$' '$\Delta t = 0.1$'},'Interpreter','latex','FontSize',14 )

        
        subplot(2,1,2)
        plot(t_sim, x_d_interp, '.','MarkerSize',2);
        title('$\Delta t = 0.01$','Interpreter','latex','FontSize',14)
        
        figure('Position',[700 500 600 400])
        subplot(2,1,1)
        plot(t_plan(1:end-1), u_d, '.');
        title( {'$u_d$' '$\Delta t = 0.1$'},'Interpreter','latex','FontSize',14 )
        
        subplot(2,1,2)
        plot(t_sim(1:end-1), u_d_interp, '.','MarkerSize',2);
        title('$\Delta t = 0.01$','Interpreter','latex','FontSize',14)
    end
end