function [T_0i, w, v] = get_T_w_v(q, qd, DH, n)

z0 = [0; 0; 1];

w = cell(n, 1); 
T_0i = cell(n, 1);
v = cell(n+1, 1);
q = [q; 0];

% Forward recursion
for i = 1:n
    R = R_DH(DH(i,:), q(i));
    p_prev = p_DH(DH(i,:));
    p = p_DH(DH(i+1,:)); % i+1 ??
    T = [R p_prev; 0 0 0 1];

    bi = z0;

    if i == 1
        T_0i{i} = T;
        w{i} = bi * qd(i);
        v{i} = [0; 0; 0];
        %alpha{i} = bi * qdd(i);
        %a{i} = R' * (-g) + cross(alpha{i}, r_i_ci{i}) + cross(w{i}, cross(w{i}, r_i_ci{i}));
        %a_e{i} = R' * (-g) + cross(alpha{i}, p) + cross(w{i}, cross(w{i}, p));
    else
        T_0i{i} = T_0i{i-1} * T;
        w{i} = R' * w{i-1} + bi * qd(i);
        v{i} = v{i-1} + cross(w{i}, p); 
        %alpha{i} = R' * alpha{i-1} + bi * qdd(i) + R' * cross(w{i-1}, bi * qd(i));
        %a{i} = R' * a_e{i-1} + cross(alpha{i}, r_i_ci{i}) + cross(w{i}, cross(w{i}, r_i_ci{i}));
        %a_e{i} = R' * a_e{i-1} + cross(alpha{i}, p) + cross(w{i}, cross(w{i}, p));
    end 
end

function R = R_DH(dh, q)
    c1 = cos(dh(1)); c2 = cos(dh(4) + q); 
    s1 = sin(dh(1)); s2 = sin(dh(4) + q);
    
    R = [c2    -s2     0 ; 
         c1*s2 c1*c2 -s1 ; 
         s1*s2 c2*s1  c1 ]; 
end
function p = p_DH(dh)    
    p = [dh(2); -sin(dh(1)) * dh(3); cos(dh(1)) * dh(3)];
end

end