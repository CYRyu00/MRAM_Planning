function object = get_u_shape_pipe()
  V = get_vertices_of_pipe([1,1,1]);
  sp1 = SupportFunction(V, 16);
  R=[0,-1,0;1,0,0;0,0,1];
  sp2 = SupportFunction(R*V, 16);
  convex_shapes = {sp1, sp2, sp1};
  c1 = [ 0.5;     1/6;0];
  c2 = [ 0.0;-0.5+1/6;0];
  c3 = [-0.5;     1/6;0];
  centers = [c1,c2,c3];
  object = ObjectGeometry(convex_shapes, centers);
end

