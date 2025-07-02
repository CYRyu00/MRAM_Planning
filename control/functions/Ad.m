function out = Ad(R, p)
out = [R zeros(3,3); S(p)*R R];
end