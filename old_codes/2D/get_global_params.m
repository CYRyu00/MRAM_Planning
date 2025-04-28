function  [m1, m2, lp, lg, m0, I0,mu,r,d,g,c_cart,c_pole,thrust_limit] = get_global_params()
global params
m1 = params(1);
m2 = params(2);
lp = params(3);
lg = params(4);
m0 = params(5);
I0 = params(6);
mu = params(7);
r = params(8);
d = params(9);
g = params(10);
c_cart = params(11);
c_pole = params(12);
thrust_limit = params(13);