function define_global_params()
% Parameters
I = [16.571710 0.830806 0.718277;
     0.830806 16.655602 1.800197;
     0.718277 1.800197 29.261652]*10e-6;
I0 = I(2,2);
thrust_limit = 0.15 ;
m0 = 0.028;
mu = 0.005964552;
r = 0.092*sqrt(2)/4;%0.032
d = 0.1;
lg = 0.08; % custom
g=9.81;
m1 = 0.2; m2 = 0.1; lp = 0.2; l1 = lp+lg;
c_cart = 100e-2; 
c_pole = 100e-5; 
global params
params = [m1, m2, lp, lg, m0, I0,mu,r,d,g,c_cart,c_pole,thrust_limit];