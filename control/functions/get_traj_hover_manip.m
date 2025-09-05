function [X_des, Xd_des, Xdd_des, Xddd_des, R_e_des, w_e_des, wd_e_des] = get_traj_hover_manip(X_hover, rpy_hover, velocity, maximum_X, N_sim, dt_sim)
    X_des = []; Xd_des = []; Xdd_des = []; Xddd_des = []; 
    R_e_des = cell(N_sim, 1); w_e_des = []; wd_e_des = [];

    X_prev = X_hover; Xd_prev = zeros(3,1); Xdd_prev = zeros(3,1);  
    R_e_prev = rpy2rot(rpy_hover(1), rpy_hover(2), rpy_hover(3));
    velocity_input = velocity;
    cnt_x = 0;
    for i = 1:N_sim
        X = X_prev + velocity * dt_sim;
        if abs(X(1)) > maximum_X(1)
            cnt_x = cnt_x + 1;
            velocity(1) = velocity_input(1) * (0.5 * cos( min(pi, cnt_x * dt_sim/0.5)) + 0.5);
            %X(1) = sign(X(1)) * maximum_X(1);
            %Xd(1) = 0;
        end    
        if abs(X(2)) > maximum_X(2)
            X(2) = sign(X(2)) * maximum_X(2);
            Xd(2) = 0;
        end
        if abs(X(3)) > maximum_X(3)
            X(3) = sign(X(3)) *maximum_X(3);
            Xd(3) = 0;
        end
        Xd = velocity;
        
        if i == 1
            Xdd = [0; 0; 0];
        else
            Xdd = (Xd - Xd_prev) / dt_sim; 
        end
        if i < 2
            Xddd = [0; 0; 0];
        else
            Xddd = (Xdd - Xdd_prev) / dt_sim;
            Xdd_prev = Xdd; 
        end

        X_des = [X_des, X];
        Xd_des = [Xd_des, Xd];
        Xdd_des = [Xdd_des, Xdd];
        Xddd_des = [Xddd_des, [0;0;0]];
        X_prev = X; Xd_prev = Xd;

        R_e = rpy2rot(rpy_hover(1), rpy_hover(2), rpy_hover(3));
        Rd = (R_e - R_e_prev) / dt_sim;
        w_e = vee(R_e' * Rd);
        R_e_prev = R_e;
        R_e_des{i} = R_e;
        w_e_des = [w_e_des, w_e];
        wd_e_des = [wd_e_des, [0;0;0]];
    end
end