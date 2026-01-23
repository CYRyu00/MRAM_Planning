mass_obj = 10;
mu_s = 0.6;
mu_d = 0.4; 

v_s = 0.2;
v_arr = -1:0.01:1;
F_n = mass_obj * 9.81;

for i = 1:length(v_arr)
    F_arr(i) = F_fric(v_arr(i), F_n, v_s, mu_s, mu_d);
end

figure;
plot(v_arr , F_arr)
grid on
function out = F_fric(v, F_n, v_s, mu_s, mu_d)
    out = F_n * (mu_d * tanh(4*v / v_s) ...
            + (mu_s - mu_d)*v/v_s/ ((v/2/v_s)^2 + 0.75)^2 );
end