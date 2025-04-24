function [ A_arr, B_arr ] = check_ctrb(x_d_interp, u_d_interp, N_sim, A_func, B_func, delta_inertia, delta_k, disturb)
    fprintf("Check Controllability for all x_d, u_d \n")
    A_arr = cell(1,N_sim);
    B_arr = cell(1,N_sim);

    for i=1:N_sim
        A_val = full(A_func(x_d_interp(i,:), u_d_interp(i,:), delta_inertia, delta_k, disturb ));  % nx x nx
        B_val = full(B_func(x_d_interp(i,:), u_d_interp(i,:), delta_inertia, delta_k, disturb ));  % nx x nu
        A_arr{i} = A_val; B_arr{i} = B_val;
        C = ctrb(A_val, B_val);
        rank_C = rank(C);
        is_controllable = (rank_C == size(A_val,1));
        if is_controllable == false
            break
        end
    end
    if is_controllable
        disp('✅ System is controllable.');
    else
        fprintf('❌ System is NOT controllable. at time step : %d\n\n',i)
    end
end

