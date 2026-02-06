classdef ObjectGeometry < handle
  properties
    support_functions
    centers
    size
  end
  
  methods
    function obj = ObjectGeometry(support_funtions, centers)
      obj.support_functions = support_funtions;
      obj.centers = centers;
      obj.size = size(support_funtions, 2);
    end    
    function s = surface(obj, idx, x)
      s = obj.support_functions{idx}.surface(x) + obj.centers(:, idx);
    end
    function s = dsurface(obj, idx, x)
      s = obj.support_functions{idx}.dsurface(x);
    end
  end
end

