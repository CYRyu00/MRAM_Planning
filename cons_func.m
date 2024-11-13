function [c,ceq]=cons_func(X_opt)
    
    global t_f
    global x_f

    x_0 = [0,0,0,0];
    
    N_st = 4;N_u=4;
    N = (length(X_opt))/(N_st+N_u);
    x = X_opt(1:N*N_st);
    u = X_opt(N*N_st+1:end);
    x = reshape(x,N,N_st);%row vectors
    u = reshape(u,N,N_u);%row vectors
    %t = linspace(0,1,N)'*t_f;
    dt = t_f/(N-1);
    
    %set initial condition and final condition
    ceq = [ x(1,:)'-x_0'; x(N,:)' - x_f']; 
    c =[];
    for k = 1:N-1
        %column vectors
        x_k = x(k,:)';
        x_k1 =  x(k+1,:)';
        x_c = (x_k + x_k1)/2;
        
        u_k = u(k,:)';
        u_k1 =  u(k+1,:)';
        u_c = (u_k + u_k1)/2;
        
        delta_k = (x_k - x_k1) + dt/6* [ FDynamics(x_k,u_k) ... 
            + 4*FDynamics(x_c,u_c) + FDynamics(x_k1,u_k1)];

        ceq = [ceq ; delta_k];

        %ineq = [u_k(5)-u_k1(5) - theta_dot_bound*dt; - u_k(5) + u_k1(5) - theta_dot_bound*dt];
        %c = [c;ineq];
    end    
    
end


