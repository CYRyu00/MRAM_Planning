function visualize_3d_cone()
    % VISUALIZE_3D_CONE
    % 3차원 공간에서 임의의 축을 기준으로 Cone 샘플링을 시각화합니다.

    clc; clear; close all;

    % --- 1. 사용자 설정 ---
    target_axis = [1; 0; 0];  % [1, 1, 1] 방향 (대각선)
    cone_angle = 10;          % 반경 30도
    num_samples = 1000;       % 샘플 개수
    
    % --- 2. 샘플링 함수 실행 ---
    % (아래에 정의된 함수 호출)
    samples = sample_3d_cone_process(target_axis, cone_angle, num_samples);

    % --- 3. 시각화 (Visualization) ---
    figure('Color', 'w', 'Position', [100, 100, 800, 600]);
    hold on; grid on; axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(135, 30); % 보기 좋은 각도로 카메라 설정
    title(sprintf('3D Cone Sampling (Axis: [%.1f, %.1f, %.1f], Angle: %d^o)', ...
        target_axis(1), target_axis(2), target_axis(3), cone_angle));

    % (A) 단위 구 (참조용 투명 구) 그리기
    [sx, sy, sz] = sphere(50);
    surf(sx, sy, sz, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % (B) 목표 중심축 그리기 (굵은 초록색 화살표)
    c_norm = target_axis / norm(target_axis);
    quiver3(0, 0, 0, c_norm(1)*1.2, c_norm(2)*1.2, c_norm(3)*1.2, ...
        'g', 'LineWidth', 3, 'MaxHeadSize', 0.5, 'DisplayName', 'Target Axis');

    % (C) 샘플링된 점 그리기 (파란색 점)
    % 경계선(Boundary)과 내부(Interior)가 섞여 있음
    scatter3(samples(1,:), samples(2,:), samples(3,:), 15, 'b', 'filled', ...
        'MarkerFaceAlpha', 0.6, 'DisplayName', 'Samples');

    % (D) 원점 표시
    plot3(0, 0, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    
    legend('show', 'Location', 'best');
    hold off;
end

function samples = sample_3d_cone_process(center_axis, angle_deg, num_samples)
    % 3D 전용 샘플링 함수 (Householder Rotation 사용)
    
    d = 3; % 3차원
    alpha_rad = deg2rad(angle_deg);
    c = center_axis(:) / norm(center_axis); % 목표 축 정규화

    % ------------------------------------------------
    % [Step 1] North Pole (Z축, [0;0;1]) 기준 로컬 샘플링
    % ------------------------------------------------
    % 3차원에서는 "Archimedes' Hat-Box Theorem"에 의해
    % Z값(cos theta)을 균일하게 뽑으면 구면 위에서도 균일해집니다.
    
    % 경계(Boundary)에 집중할 비율 (예: 80%)
    n_bound = floor(num_samples * 0.8);
    n_inside = num_samples - n_bound;
    
    % (A) 내부 샘플: cos(theta)를 [cos(alpha), 1] 사이에서 균일하게 뽑음
    z_min = cos(alpha_rad);
    z_inside = z_min + (1 - z_min) * rand(1, n_inside);
    
    % (B) 경계 샘플: cos(theta)를 cos(alpha)로 고정
    z_bound = ones(1, n_bound) * z_min;
    
    % Z값 합치기
    z = [z_bound, z_inside];
    
    % X, Y 좌표 계산 (반지름 r = sqrt(1 - z^2))
    r = sqrt(1 - z.^2);
    phi = 2 * pi * rand(1, num_samples); % 방위각은 0~360도 균일
    
    x = r .* cos(phi);
    y = r .* sin(phi);
    
    % 로컬 샘플 조립 (기본 축은 Z축 [0;0;1] 이라고 가정)
    % u_local = [x; y; z]
    u_local = [x; y; z];

    % ------------------------------------------------
    % [Step 2] 회전 (Rotation): Z축([0;0;1]) -> 목표 축(c)
    % Householder Reflection 사용
    % ------------------------------------------------
    e3 = [0; 0; 1]; % 로컬 기준 축
    
    % (예외 처리) 이미 축이 같으면 회전 안 함
    if norm(c - e3) < 1e-8
        samples = u_local;
        return;
    end
    if norm(c + e3) < 1e-8
        samples = -u_local;
        return;
    end
    
    % 반사 벡터 h = e3 - c
    h = e3 - c;
    h_sq = h' * h;
    
    % H * x = x - 2 * h * (h' * x) / (h' * h)
    dot_vals = h' * u_local; % (1 x N)
    samples = u_local - (2 / h_sq) * (h * dot_vals);
end