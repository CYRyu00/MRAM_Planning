function [Sa] = S(a)
   Sa = [0, -a(3), a(2);
         a(3), 0, -a(1);
         -a(2), a(1), 0];
end

