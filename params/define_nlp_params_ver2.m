function [Q1, Q2, Q3, R] = define_nlp_params_ver2(nu)
Q1 = diag([0, 0, 0, 0, 1, 1, 1, 1]) * 1e-1; % q_dot
Q2 = diag([1, 1, 1, 0]) * 1e0; % q_o_ref 
Q3 = diag([1, 1, 1, 1, 1, 1, 1, 1]) * 1e1;% only for hover
R = diag(ones(nu, 1)) * 1e0; % u
end
