clc; clear; close all;
addpath("../dynamics", "../functions", "../../params" )
clear; close all
params = define_params_ver2(); 
mu = params{3}; r = params{4};
mb = params{17}; cb = params{18}; Ib = params{19};
mt = params{20}; ct = params{21}; It = params{22};
ma = params{23}; ca = params{24}; Ia = params{25};

% totoal robotic arm
Ib = Ia + Ib;% Ib = Ib * 0.01;
mb = ma + mb;
m0 = mt + mb;

B_tau = [-r r];
B_f = [1 1];

e_1 = [1; 0; 0];
e_2 = [0; 1; 0];
e_3 = [0; 0; 1];
%% inertia
num_AMs = 2;
AM_mass = 0; % mass of shape
AM_inertia = zeros(3, 3); % inertia w.r.t. its com
AM_com = [0; 0; 0];% 0 to com
r_0j = cell(num_AMs, 1); % 0 to j'th module
r_cj = cell(num_AMs, 1); % com to j'th module
I_cj = cell(num_AMs, 1); % inertia of j'th module w.r.t. com of shpe
mass_ams = m0 * ones(num_AMs, 1);
R_shape = cell(1, num_AMs);
% shape_pitch =[0 -10 -20 -10 0 0 10];
% shape_pitch =[0 -10 -10 10 -10 -10 10];
shape_pitch =[0 0 0 0 0 0 0];

for j = 1:length(shape_pitch)
    R_shape{j} = Ry(shape_pitch(j) / 180 *pi);
end
l1 = 0.35; l2 = 0.30; % ca = 0.24

AM_mass = sum(mass_ams);
for j = 1:num_AMs
    if j == 1
        r_0j{j} = [0; 0; 0];
    else
        r_0j{j} = r_0j{j-1} -l2 * R_shape{j-1} * e_1 - l1 * R_shape{j} * e_1;% p_core to j
    end
    AM_com = AM_com + r_0j{j} * mass_ams(j)/AM_mass;
end
for j = 1:num_AMs
    r_cj{j} = r_0j{j} - AM_com;% p_com to j
end

%% 1. 시스템 정의 (B 행렬 설정)
% w = [fx; fz; tau_y]

B = zeros(3, num_AMs*2);
theta_arr = [10, -10, 10, -10, 20, -20];
for j = 1:num_AMs
    uj = [sin(theta_arr(j) /180*pi); 0; cos(theta_arr(j) /180*pi)]; % TODO
    B(:, 2*j-1:2*j) = [[e_1, e_3]' * uj * B_f; e_2'*S(r_cj{j}) * uj * B_f + B_tau];
end

n_actuators = 2*num_AMs; % 입력 개수
f_min = 0;
f_max = 20;

%(Vertex Enumeration)
num_vertices = 2^n_actuators;
perms = dec2bin(0:num_vertices-1) - '0';
f_vertices = perms' * (f_max - f_min) + f_min;

%Wrench 공간으로 매핑 (Linear Mapping)
w_points = B * f_vertices; 
% w_points = w_points - m0* 9.81 *e_2;

% 행 1: f_x, 행 2: f_z, 행 3: tau_y

% 볼록 껍질 (Convex Hull) 계산
x = w_points(1, :)';
z = w_points(2, :)';
y = w_points(3, :)';

[K, volume] = convhull(x, y, z);

% 내접구 반지름 (Chebyshev Radius) 계산
% 원점(0,0,0)에서 각 면(Face)까지의 수직 거리를 계산하여 최소값을 찾습니다.
% 이것이 곧 '어느 방향으로든 낼 수 있는 최소한의 힘/토크'를 의미합니다.

min_dist = inf;

for i = 1:size(K, 1)
    % 삼각형 면을 이루는 세 점 추출
    p1 = w_points(:, K(i, 1));
    p2 = w_points(:, K(i, 2));
    p3 = w_points(:, K(i, 3));
    
    % 삼각형의 법선 벡터 (Normal Vector) 계산
    normal_vec = cross(p2 - p1, p3 - p1);
    normal_vec = normal_vec / norm(normal_vec); % 단위 벡터화
    
    % 원점에서 평면까지의 거리 공식: d = |n . p1|
    % (평면 방정식이 n.(x - p1) = 0 이므로)
    dist = abs(dot(normal_vec, p1));
    
    if dist < min_dist
        min_dist = dist;
    end
end

%% 결과 출력 및 시각화
fprintf('=== Wrench Polytope Analysis ===\n');
fprintf('Polytope Volume (기동성 총량): %.4f\n', volume);
fprintf('Inscribed Sphere Radius (강건성/여유분): %.4f\n', min_dist);

figure('Color', 'w');
% 볼록 껍질 그리기 (투명도 alpha 적용)
trisurf(K, x, y, z, 'FaceColor', 'cyan', 'FaceAlpha', 0.5, 'EdgeColor', 'b');
hold on;

% 원본 꼭짓점 점 찍기
plot3(x, y, z, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');

% 좌표축 및 라벨 설정
grid on; axis equal;
xlabel('Force X (f_x)');
zlabel('Force z (f_z)');
ylabel('Torque (\tau)');
title(['Wrench Polytope (Vol: ', num2str(volume, '%.2f'), ...
       ', Radius: ', num2str(min_dist, '%.2f'), ')']);

% 원점 표시
plot3(0,0,0, 'k+', 'MarkerSize', 15, 'LineWidth', 2);
% view(1); % 3D 뷰 설정

