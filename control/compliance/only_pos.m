clear; close all;
%% Dynamics Parameters
n_am = 4;
mass = 2.0 * ones(1,n_am); m_com = sum(mass);
k_spring = 300; k_s = [0, ones(1, n_am - 1) 0] * k_spring; % N/m , spring coeff
l1 = 0.3; l2 = 0.2;
pc_i = zeros(3, n_am);
for i = 1:n_am
    pc_i(1, i) = -(l1 + l2) * (i - (n_am+1)/2);
end
e_1 = [1;0;0]; e_2 = [0;1;0]; e_3 = [0;0;1];
g = 9.81;

%% Simulation parameters
dt = 0.001; N = 10/ dt; 
show_video = true;
save_video = false;
video_speed = 2;
gap = (l1+l2) *2;

%% Control parameters
wn = 3; damp = 1.1; % 0.5, 1.2
kp = wn^2; 
kv = 2 * damp *sqrt(kp);
%kp = 12; kv = 10;
alpha = 4;
gamma = alpha * 1;
beta = 0.2;
epsilon = 1;

alpha_cen = 4 * 1;
gamma_cen = alpha_cen * 0.25;
%% Uncertainty 
spring_error = 1.0;
delta = [1 0 0 0; 0 0 0 0; -1 0 0 0] * diag(mass) * g * 1;

%% Trajectory
p_com_ref = [0.1; 0 ; 0.2];
p_tilde_ref = [2 1 -1 -2; 0 0 0 0; -3 -1 1 3] * 0.01;
pd_tilde_ref = [0 0 0 0; 0 0 0 0; -2 -1 1 2] * 0.0;
R_c_ref = eye(3);
pd_com_ref = [1; 0; 3]*0.1;
w_c_ref = [0; 0.00; 0];
[p_des_traj, pd_des_traj, pdd_des_traj, pddd_des_traj, P_des_traj, Pd_des_traj, Pdd_des_traj, Pddd_des_traj ] ...
           = gen_traj_hover(p_com_ref, pd_com_ref, p_tilde_ref, pd_tilde_ref, R_c_ref, w_c_ref , pc_i, n_am, N, dt);
%% Backstepping Controller
M = diag(mass); M = kron(M, eye(3));
M_bar = diag([sum(mass), mass(1:end-1)]); M_bar = kron(M_bar, eye(3));
K_bar = k_spring *[0 0 0 0; 0 1 -1 0;0 -1 2 -1;0 1 0 3]; K_bar = M_bar \ kron(K_bar, eye(3)); % todo
B_bar = [1 1 1 1; 3/4 -1/4 -1/4 -1/4; -1/4 3/4 -1/4 -1/4;-1/4 -1/4 3/4 -1/4]; B_bar = kron(B_bar, eye(3));

p_i = pc_i + ones(3, n_am) * 0.0; % p1, p2 ~
pd_i = ones(3, n_am) * 0; pdd_i = ones(3, n_am) * 0; 
R_c = eye(3,3); Rd_c = zeros(3, 3);
u = zeros(3, n_am); ud = zeros(3, n_am); % lamda R e3 - m g e3
u_cen = zeros(3 * n_am, 1); ud_cen = zeros(3 * n_am, 1);
force = zeros(3, n_am); forced = zeros(3, n_am);
delta_hat = zeros(3, n_am); delta_hatd = zeros(3, n_am);

e_p_arr = []; e_pd_arr = []; 
p_arr = []; p_i_arr = []; P_arr = [];
force_arr = []; delta_arr = [];
times = [];


for i = 1:N
    p_com = mean(p_i, 2);
    pd_com = mean(pd_i, 2);
    p_i_tilde = zeros(3, n_am);
    pd_i_tilde = zeros(3, n_am);
    
    P = p_com;
    Pd = pd_com;
    for j = 1:n_am
        p_i_tilde(:, j) = p_i(:, j) - p_com - R_c * pc_i(:, j);
        pd_i_tilde(:, j) = pd_i(:, j) - pd_com - Rd_c * pc_i(:, j);
        
        if j < n_am
            P = [P; p_i_tilde(:, j)];
            Pd = [Pd; pd_i_tilde(:, j)];
        end
    end
    
    % Central control
    e_p_bar = P - P_des_traj(:, i);
    e_pd_bar = Pd - Pd_des_traj(:, i);
    K_bar_hat = K_bar * spring_error;
    nu_e_cen = M_bar\B_bar * u_cen - K_bar_hat * P_des_traj(:, i) + 0; % todo rotational spring force estimation    
    
    eta_cen = -alpha_cen * nu_e_cen - gamma_cen * (e_pd_bar + epsilon * e_p_bar);
    nud_cen = K_bar_hat * zeros(3*n_am, 1) - 0; % todo Pd_des, rotational spring force estimation;
    
    ud_cen = B_bar\M_bar * (nud_cen + eta_cen); 
    u_cen = u_cen + ud_cen * dt;

    delta_s_hat = -(M*B_bar)\M_bar * K_bar_hat * P; 
    delta_s_hat = reshape(delta_s_hat, 3, n_am);

    % Decentralized control
    for j = 1:n_am
        pddd_des = pddd_des_traj(3*j-2:3*j, i);
        pdd_des = pdd_des_traj(3*j-2:3*j, i);
        pd_des = pd_des_traj(3*j-2:3*j, i);
        p_des = p_des_traj(3*j-2:3*j, i);
        
        e_p = p_i(:, j) - p_des;
        e_pd = pd_i(:, j) - pd_des;
        pdd_i_hat = force(:, j)/mass(j) - g * e_3 + delta_hat(:, j) + delta_s_hat(:, j); % todo spring
        %pdd_i_hat = pdd_i(:, j);
        e_pdd_hat = pdd_i_hat - pdd_des;

        % generation error 
        nu_e = u(:, j)/mass(j) - pdd_des + kp*e_p + kv*e_pd + delta_hat(:, j);
        % adaptive law
        delta_hatd(:, j) = beta * (e_pd + epsilon * e_p + kv /gamma * nu_e );
        
        eta = -alpha * nu_e - gamma * (e_pd + epsilon * e_p);
        nud = pddd_des - kv * e_pdd_hat - kp*e_pd - delta_hatd(:, j);
        ud(:, j) = mass(j) * (nud + eta); 
        u(:, j) = u(:, j) + ud(:, j) *dt;
        delta_hat(:, j) = delta_hat(:, j) + delta_hatd(:, j) *dt;
    end
    
    force = u + reshape(u_cen, 3, n_am) + [0;0;1] * mass * g;
    forced = ud + reshape(ud_cen, 3, n_am);

    input_FD = force + delta;

    pdd_i = ForwardDynamics(input_FD, p_i, R_c, n_am, mass, g, k_s, pc_i);
    p_i = p_i + pd_i * dt + 0.5 * pdd_i *dt *dt;
    pd_i = pd_i + pdd_i * dt;

    e_p_arr = [e_p_arr, e_p];
    e_pd_arr = [e_pd_arr, e_pd]; 
    p_arr = [p_arr, p_i(:, end)];
    p_i_arr = [p_i_arr, p_i(:)];
    P_arr = [P_arr, P];

    force_arr = [force_arr, force(:)]; 
    delta_arr = [delta_arr, delta(:)];

    times = [times, i*dt];
end
%% plot
figure
subplot(3,1,1)
plot(times, p_arr)
legend("x","y","z")
subplot(3,1,2)
plot(times, e_p_arr)
legend("x","y","z")
subplot(3,1,3)
plot(times, e_pd_arr)
legend("x","y","z")

%% Video
if show_video
figure('Position',[600 100 500 500]);
dN = 0.1 / dt;
framesPerSecond = 1/dt/dN * video_speed;
rate = rateControl(framesPerSecond);
arrow_len = 0.05;

if save_video
    video_filename = '../images/test.avi';
    video = VideoWriter(video_filename);
    video.FrameRate = framesPerSecond;
    open(video);
end

for i = 1:dN:N
    clf;
    grid on; axis equal;
    set(gca, 'GridLineWidth', 1.0, 'GridAlpha', 0.3);
    xmin = p_des_traj(end-2, i) - gap;
    xmax = p_des_traj(1, i) + gap;
    xlim([xmin , xmax])
    
    ymin = p_des_traj(3, i) - gap;
    ymax = p_des_traj(end, i) + gap;
    ylim([ymin, ymax])
    
    hold on;
    R = eye(3);
    
    for j = 1:n_am
        R_e_j = eye(3);
        X_quad = p_i_arr(3*j-2 : 3*j, i);
        R_quad = Ry(0);
        plot(X_quad(1), X_quad(3), 'o', 'Color', 'b');
        
        dx = force_arr(3*j-2, i) * mass(j) / g * arrow_len;
        dz = force_arr(3*j, i) * mass(j) / g * arrow_len;
        quiver(X_quad(1), X_quad(3), dx, dz, 0, 'r', 'LineWidth', 1.5);
        
        X1 = X_quad + l1 * R_e_j * e_1;
        X2 = X_quad - l2 * R_e_j * e_1;
        plot([X1(1), X2(1)], [X1(3), X2(3)], 'b-', 'LineWidth', 1.0);

        % external force
        dx = delta_arr(3*j-2, i) * mass(j) / g * arrow_len;
        dz = delta_arr(3*j, i) * mass(j) / g * arrow_len;
        quiver(X_quad(1), X_quad(3), dx, dz, 0, 'g', 'LineWidth', 1.5);
    end

    for j = 1:n_am
        R_e_des_j = eye(3);
        X_des_quad = p_des_traj(3*j-2:3*j, i);      
        plot(X_des_quad(1), X_des_quad(3), 'o', 'Color', 'black', 'LineWidth', 2.0);

        X1 = X_des_quad + l1 * R_e_des_j * e_1;
        X2 = X_des_quad - l2 * R_e_des_j * e_1;
        plot([X1(1), X2(1)], [X1(3), X2(3)], 'k--', 'LineWidth', 2.0);
    end

    % com
    plot(P_arr(1, i), P_arr(3, i), 'o', 'Color', 'b', 'MarkerSize', 10);
    plot(P_des_traj(1, i), P_des_traj(3, i), 'o', 'Color', 'k', 'MarkerSize', 8);

    xlabel('X position'); ylabel('Z position');
    title_string = sprintf("time : %.2f sec", times(i));
    title(title_string);
    drawnow

    if save_video
        frame = getframe(gcf);
        writeVideo(video, frame);
        waitfor(rate);
    end
end
if save_video
    close(video)
    fprintf("\n Video saved at %s\n", video_filename);
end
end
%%
function pdd_i = ForwardDynamics(u, p_i, R_c, n_am, mass, g, k_s, pc_i)
    q = zeros(3, n_am + 1); % q(3, 1) is 0 always
    p_com = mean(p_i, 2);
    p_i_tilde = zeros(3, n_am);
    for i = 1:n_am
        p_i_tilde(:, i) = p_i(:, i) - p_com - R_c * pc_i(:, i);
        if i > 1
            q(:, i) = p_i_tilde(:, i) - p_i_tilde(:, i-1);
        end
    end

    pdd_i = zeros(3, n_am);
    for i = 1:n_am
        pdd_i(:, i) = -g * [0;0;1] + (u(:, i) - k_s(i) * q(:, i) + k_s(i+1) * q(:, i+1)) / mass(i);
    end

end

function out = Rx(q)
    out = [1, 0, 0;
          0, cos(q), -sin(q);
          0, sin(q), cos(q)];
end
function out = Ry(q)
    out = [cos(q), 0, sin(q);
           0, 1, 0;
           -sin(q), 0, cos(q)];
end