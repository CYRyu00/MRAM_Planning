function [p_des_traj, pd_des_traj, pdd_des_traj, pddd_des_traj, P_des_traj, Pd_des_traj, Pdd_des_traj, Pddd_des_traj ] ...
           = gen_traj_hover(p_com_ref, pd_com_ref, p_tilde_ref, pd_tilde_ref, R_c_ref, w_c_ref , pc_i, n_am, N, dt)
    
% p1, p2 ...
p_des_traj = zeros(3*n_am, N);
pd_des_traj = zeros(3*n_am, N);
pdd_des_traj = zeros(3*n_am, N);
pddd_des_traj = zeros(3*n_am, N);

% p_com, p1 tilde, p2 tilde ...
P_des_traj = zeros(3*n_am, N);
Pd_des_traj = zeros(3*n_am, N);
Pdd_des_traj = zeros(3*n_am, N);
Pddd_des_traj = zeros(3*n_am, N);

p_com = p_com_ref;
R_c = R_c_ref;
Rd_c = R_c * S(w_c_ref);
for i = 1:N       
    P_des_traj(1:3, i) = p_com;
    Pd_des_traj(1:3, i) = pd_com_ref;

    for j = 1:n_am-1
        P_des_traj (3*j+1:3*j+3, i) = p_tilde_ref(:, j) + pd_tilde_ref(:, j) * (i - 1) * dt;
        Pd_des_traj (3*j+1:3*j+3, i) = pd_tilde_ref(:, j);
    end
    for j = 1:n_am
        p_des_traj (3*j-2:3*j, i) = p_com + R_c * pc_i(:, j) + p_tilde_ref(:, j) + pd_tilde_ref(:, j) * (i - 1) * dt;
        pd_des_traj (3*j-2:3*j, i) = pd_com_ref + Rd_c * pc_i(:, j) + pd_tilde_ref(:, j);
    end
    
    p_com = p_com + pd_com_ref *dt;
    Rd_c = R_c * S(w_c_ref);
    R_c = R_c + Rd_c * dt;
end
end

function [Sa] = S(a)
   Sa = [0, -a(3), a(2);
         a(3), 0, -a(1);
         -a(2), a(1), 0];
end


