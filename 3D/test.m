%%
m = 10;
tic
shapes = generate_all_shapes(m);
toc
%%
shape = [...
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     1     1     1     1     0     0     0     0     0     0;
     2     1     1     1     0     0     0     0     0     0; %center
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0];

%shape = shapes{10}{2};
for i=1:length(shapes)
    for j=1:length(shapes{i})
        if shape == shapes{i}{j}
           disp([i,j])
        end
    end
end

params = define_params();
m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5};
[AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);


%%
num_shapes = zeros(m);
for i=1:length(shapes)
    num_shapes(i) = length(shapes{i});
    num_shapes(i) = num_shapes(i) / 3.7^(i-1);
end
plot(num_shapes)

%% generate door Example Usage:
n=4;
dh = [0,0,0.95,0;   % [alpha, a, d, theta]
      -pi/2, 0.9 , 0,0;
      0,-0.1,0.23,pi/2;
      pi/2,0,0,-pi/2;
      pi/2,0,0,0];
g = [0;0;-9.81];

params = define_params();
m0 = params{1}; I0 = params{2};mu = params{3}; r= params{4}; d= params{5};thrust_limit= params{6};

m=7;
all_shapes = generate_all_shapes(m);
shape = all_shapes{7}{10};
[AM_com, AM_mass, AM_inertia] = get_inertia(shape ,m0, I0, d);
mass = {20, 3, 0.5, AM_mass , 0};
inertia = {eye(3)*1, eye(3)*2, eye(3)*3, AM_inertia, zeros(3,3)};
r_i_ci = {[0.5;-0.02;0.05],[-0.05;0;0.08],[0;0;-0.05],[AM_com(1);0;AM_com(2)], zeros(3,1)};

do_view =1; q=[pi/6;-pi/3;pi/3;-pi/4];
robot = generate_door_ver2(n,dh,r_i_ci, g, mass,inertia, do_view,q);