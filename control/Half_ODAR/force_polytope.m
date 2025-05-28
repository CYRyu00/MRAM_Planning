addpath("../dynamics", "../functions", "../../params" )
clear; close all
%%
params = define_params();
m0 = params{1}; I0 = params{2}; mu = params{3}; r= params{4}; d= params{5};
thrust_limit= params{6}; kt = params{7}; c_1 = params{8}; c_2 = params{9}; mass_door = params{10};
handle_factor = params{11}; inertia = params{12}; r_i_ci = params{13}; n = params{14}; dh = params{15}; gravity = params{16};
    
eta_arr = ([1, 5:5:55]) / 180 * pi;
num_figs = length(eta_arr);
cnt = 0;

% 화면 정렬 설정
cols = 4;  % 한 줄에 최대 4개
rows = ceil(num_figs / cols);

fig_width = 400;   
fig_height = 250;  
x_gap = 20;        
y_gap = 90;      
screen_offset_x = 50;
screen_offset_y = 50;

for eta = eta_arr
    cnt = cnt + 1;
    s = sin(eta); c = cos(eta);
    A_eta = [mu*s/sqrt(2) - r*c, mu*s/sqrt(2) - r*c, - mu*s/sqrt(2) + r*c, - mu*s/sqrt(2) + r*c;
             - mu*s/sqrt(2) + r*c, mu*s/sqrt(2) - r*c, mu*s/sqrt(2) - r*c, - mu*s/sqrt(2) + r*c;
             mu*c + sqrt(2)*r*s, - mu*c - sqrt(2)*r*s, mu*c + sqrt(2)*r*s, - mu*c - sqrt(2)*r*s;
             s/sqrt(2), - s/sqrt(2), - s/sqrt(2), s/sqrt(2);
             - s/sqrt(2), - s/sqrt(2), s/sqrt(2), s/sqrt(2);
             c, c, c, c];
    fprintf("eta : %.1f\n", eta /pi *180)
    fprintf("rank A_eta: %d\n", rank(A_eta));
    disp(A_eta)
    D = 2*d;
    p=[];
    p(:,1) = [0; 0; 0]; p(:,2) = [D; 0; 0]; % p(:,3) = [D; 0; 0];
    p_avg = mean(p, 2);
    R{1} = eye(3,3); R{2} = [0 -1 0;1 0 0; 0 0 1]; % R{3} = eye(3,3);
    
    A_shape = [];
    for i=1:length(p(1,:))
        A_shape = [A_shape, Ad(R{i}, p(:,i) - p_avg) * A_eta];
    end
    fprintf("rank A_shape: %d\n\n", rank(A_shape));
    disp(A_shape)
    
    A_tau = A_shape(1:3, :);
    A_force = A_shape(4:6, :);
    
    thrusts = dec2bin(0:255) - '0';  % size: 256 x 8
    thrusts = 2 * thrusts - 1;     
    
    valid_forces = []; 
    
    for i = 1:size(thrusts, 1)
        u = thrusts(i, :)';
        tau = A_tau * u;
        if norm(tau) < 1e-6 
            f = A_force * u;
            valid_forces = [valid_forces; f'];
        end
    end
    
    % Convex hull of valid forces
    if isempty(valid_forces)
        disp('No valid force found that satisfies tau = 0');
    else
        K = convhull(valid_forces(:,1), valid_forces(:,2), valid_forces(:,3));

        col_idx = mod(cnt-1, cols);
        row_idx = floor((cnt-1) / cols);
    
        pos_x = screen_offset_x + col_idx * (fig_width + x_gap);
        pos_y = screen_offset_y + (rows - row_idx - 1) * (fig_height + y_gap);  % y는 위에서 아래로
    
        figure(cnt)
        set(gcf, 'Position', [pos_x, pos_y, fig_width, fig_height]);
        axis equal
        trisurf(K, valid_forces(:,1), valid_forces(:,2), valid_forces(:,3), ...
            'FaceAlpha', 0.6, 'EdgeColor', 'k');
        xlabel('$F_x$', 'Interpreter', 'latex');
        ylabel('$F_y$', 'Interpreter', 'latex');
        zlabel('$F_z$', 'Interpreter', 'latex');
        
        title_eta = sprintf('Force Polytope $\\tau = 0$, $\\eta = %.2f^\\circ$', eta * 180 / pi);
        title(title_eta, 'Interpreter', 'latex');

        grid on
        xlim([-3, 3]); ylim([-3, 3]); zlim([-10, 10]);
    end
end
%%
function out = Ad(R, p)
out = [R zeros(3,3); S(p)*R R];
end
function out = S(p)
out = [0 -p(3) p(2);
       p(3) 0 -p(2);
       -p(2) p(1) 0];
end