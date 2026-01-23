% 1. 시스템 정의 (예: s / (s + k), k=2)
kv = 6;
kp = 8;
ki = 10;
sys = tf([0 1 1], [1 kv kp]);
% sys = tf([1 0 0 0], [1 kv kp ki]);


fprintf('--- 시스템 분석 시작 ---\n');

%% 2. L-2 Gain (H-inf norm) 계산
l2_gain = norm(sys, inf);
fprintf('1. L-2 Gain (H-inf norm): %.4f\n', l2_gain);

%% 3. L-inf Gain (L-1 Norm) 계산
% 상태 공간 행렬 추출 (A, B, C, D)
[A, B, C, D] = ssdata(sys);

if D ~= 0
    % Case: Proper but not strictly proper (예: s/(s+k))
    % Feedthrough 항(D)이 델타 함수를 생성함
    
    % 1) Strictly proper part (D를 제외한 부분)의 임펄스 응답 계산
    sys_p = ss(A, B, C, 0);
    [h_p, t_p] = impulse(sys_p);
    
    % 2) L-inf gain = |D| + integral(|h_p(t)|)
    l_inf_gain = abs(D) + trapz(t_p, abs(h_p));
    
    fprintf('2. L-inf Gain (L-1 norm): %.4f (Feedthrough D=%.1f 포함)\n', l_inf_gain, D);
else
    % Case: Strictly proper (예: 1/(s+k))
    [h, t] = impulse(sys);
    l_inf_gain = trapz(t, abs(h));
    
    fprintf('2. L-inf Gain (L-1 norm): %.4f\n', l_inf_gain);
end