classdef ObjectGeometry < handle
  properties
    convex_shapes
    centers
    size
  end
  
  methods
    function obj = ObjectGeometry(convex_shapes, centers)
      obj.convex_shapes = convex_shapes;
      obj.centers = centers;
      obj.size = size(convex_shapes, 2);
    end    
    function s = surface(obj, idx, x)
      s = obj.convex_shapes{idx}.surface(x) + obj.centers(:, idx);
    end
    function s = dsurface(obj, idx, x)
      s = obj.convex_shapes{idx}.dsurface(x);
    end
  end
end

