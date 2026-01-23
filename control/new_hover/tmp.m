% 1. 전달 함수 정의 (예: G(s) = (s+2) / (s^2 + 2s + 2))
num = [1 0 0];
den = [1 6 8];
sys = tf(num, den);

% 2. 안정성 확인 (모든 Pole의 실수부가 음수인지)
poles = pole(sys);
is_stable = all(real(poles) < 0);

% 3. Relative Degree 확인 (0 또는 1이어야 함)
rel_deg = length(den) - length(num);
is_rel_deg_ok = (rel_deg == 0 || rel_deg == 1);

% 4. 실수부(Real Part) 양수 여부 확인
% 주파수 응답을 계산하여 최소값이 0보다 큰지 확인합니다.
[re, im, w] = nyquist(sys);
re_squeezed = squeeze(re);
min_real = min(re_squeezed);
is_real_pos = min_real > 0;

% 최종 결과 출력
if is_stable && is_rel_deg_ok && is_real_pos
    disp('이 시스템은 SPR입니다.');
else
    disp('이 시스템은 SPR이 아닙니다.');
end
nyquist(sys);
grid on;