function samples = sample_3d_cone_process(center_axis, angle_deg, num_samples, bound_ratio)
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
    n_bound = floor(num_samples * bound_ratio);
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