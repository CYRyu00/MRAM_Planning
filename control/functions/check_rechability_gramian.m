function [min_eigval_arr, max_eigval_arr,xT_W_r_x_arr ] = check_rechability_gramian(A_arr, B_arr, N_horizon, dt_sim, N_sim, n, do_print)
    min_eigval_arr = [];
    max_eigval_arr = [];
    xT_W_r_x_arr = [];
    
    for i = 1:N_sim - N_horizon + 1
        W_r = zeros(2*n);
        Phi = eye(2*n);
    
        for j = 0:N_horizon-1
            A_j = A_arr{i + j};
            B_j = B_arr{i + j};
    
            if j > 0
                Phi = A_arr{i + j - 1} * Phi;
            end
    
            W_r = W_r + Phi * B_j * B_j' * Phi';
        end
        [V, D] = eig(W_r);  % D: diagonal matrix of eigenvalues, V: eigenvectors
        [min_eigval, min_idx] = min(diag(D)); 
        [max_eigval, max_idx] = max(diag(D)); 
        min_eigval_arr = [min_eigval_arr; min_eigval];
        max_eigval_arr = [max_eigval_arr; max_eigval];
        
        if do_print
            fprintf('Step %d : \n', i);
            fprintf('Eigenvalues of W_r:\n');
            disp(diag(D)');
        
            [max_eigval, max_idx] = max(diag(D));
            fprintf('Largest eigenvalue: %g\n', max_eigval);
            fprintf('Corresponding eigenvector:\n');
            disp(V(:, max_idx)');
        
            fprintf('Smallest eigenvalue: %g\n', min_eigval);
            fprintf('Corresponding eigenvector:\n');
            disp(V(:, min_idx)');
            
            fprintf('-----------------------------------\n');
        end
        xT_W_r_x = zeros(1,2*n);
        for k = 1:2*n
            e_k = zeros(2*n,1); e_k(k) = 1;
            xT_W_r_x(k) = e_k'*W_r*e_k;
        end
        xT_W_r_x_arr = [xT_W_r_x_arr; xT_W_r_x];
    end
    
    t_tmp = (1:N_sim - N_horizon + 1)*dt_sim;
    figure('Position',[100 700 600 300])
    hold on
    plot( t_tmp, log(max_eigval_arr)/log(10))
    plot( t_tmp, log(min_eigval_arr)/log(10))
    legend("max eigen value", "min eigen value",'FontSize', 14,'Interpreter', 'latex')
    xlabel('time[sec]' , 'FontSize', 14, 'Interpreter', 'latex')
    ylabel( '$\log_{10}{(\lambda)}$','FontSize', 14,'Interpreter', 'latex')
    title({'Eigenvalues of $W_r$ (log scale)'}, ...
           'Interpreter','latex','FontSize',14);
    grid on
    
    figure('Position',[100 300 600 300])
    hold on
    colors = lines(n);
    for k=1:n
        plot( t_tmp', log(xT_W_r_x_arr(:,k)),'Color', colors(k,:))
        plot( t_tmp', log(xT_W_r_x_arr(:,k+n)),'Color', colors(k,:) , 'LineStyle','--')
    end 
    legend({'$q_1$', '$\dot{q_1}$','$q_2$','$\dot{q_2}$','$q_3$','$\dot{q_3}$','$q_4$','$\dot{q_4}$'}, ...
           'Interpreter','latex','FontSize',14);
    xlabel('time[sec]', 'FontSize', 14, 'Interpreter', 'latex')
    ylabel( '$\log_{10}{(e_i^T W_r e_i)}$','FontSize', 14,'Interpreter', 'latex')
    title({'$x^T W_r x$'}, ...
           'Interpreter','latex','FontSize',14);
    grid on
end

