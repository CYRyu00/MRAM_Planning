function [cost]= obj_func(X_opt)
    global t_f 
    global x_f

    N_st = 4;N_u=4;
    N = (length(X_opt))/(N_st+N_u);
    x = X_opt(1:N*N_st);
    u = X_opt(N*N_st+1:end);
    x = reshape(x,N,N_st);%row vectors
    u = reshape(u,N,N_u);%row vectors
    %t = linspace(0,1,N)'*t_f;
    dt = t_f/N;
    
    cost = 0;
    Q = diag([1,1,1,1])*0.1; R=diag([1,1,1,1]); % theta_tilt dot
    u = [u;u(N,:)];
    for i=1:N
        cost = cost + (x(i,:)-x_f)*Q*(x(i,:)-x_f)' + u(i,:)*R*u(i,:)';
    end    
end