classdef SupportFunction < handle 
  properties
    V
    p
    n
  end
  
  methods
    function obj = SupportFunction(V, p)
      obj.V = V;
      obj.p = p;
      obj.n = size(V, 2);
    end
    
    function s = surface(obj, x)
      dot = obj.V'*x;
      A = max(dot, zeros(obj.n, 1));
      tilde_a_p = sum(power(A, obj.p));
      hat_a_p1 = power(A, obj.p-1);
      s = obj.V * power(tilde_a_p, 1/obj.p - 1) * hat_a_p1;
    end

    function ds = dsurface(obj, x)
      dot = obj.V'*x;
      A = max(dot, zeros(obj.n, 1));
      tilde_a_p = sum(power(A, obj.p));
      hat_a_p1 = power(A, obj.p-1);
      hat_a_p2 = power(A, obj.p-2);
      ds = (obj.p-1) * obj.V * (power(tilde_a_p, 1/obj.p - 1) * diag(hat_a_p2) - power(tilde_a_p, 1/obj.p - 2) * (hat_a_p1 * hat_a_p1') ) * obj.V';
    end
  end
end

