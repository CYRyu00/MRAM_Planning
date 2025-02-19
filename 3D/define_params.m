function params = define_params()
% Parameters
if nargin < 1
        scale = 1; % Default value if scale is not provided
end

I0 = [16.571710 0.830806 0.718277;
     0.830806 16.655602 1.800197;
     0.718277 1.800197 29.261652]*10e-6;
thrust_limit = 0.15 ;
m0 = 0.028;
mu = 0.005964552;
r = 0.092*sqrt(2)/4;%0.032
d = 0.1;
lg = 0.08; % custom

c_cart = 100e-2; 
c_pole = 100e-5; 

params = {m0, I0, mu, r, d, thrust_limit, lg,c_cart,c_pole,};